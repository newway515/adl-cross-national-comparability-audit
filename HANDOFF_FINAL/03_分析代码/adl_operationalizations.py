# -*- coding: utf-8 -*-
"""
M2 操作化模块 + 全 13 人群 2PL 冒烟
========================================
采集 CHARLS 原始 4 级量表生成 M2b/M2a 评分层，
并对 13 人群逐库跑 2PL IRT 前置诊断 (information function, convergence, model fit).

CHARLS M2 操作化 (V9.1 冻结):
  - M2b (主分析): 二值化阈值口径 - "需帮助或不能做 (3 或 4) = 1", 其余=0
  - M2a (敏感性): 严重度加权 - 五条目 4 级取均值作为连续分

对其他 12 人群: M2 = M1 (无评分梯度数据, E2 不可估计)

用法:
    python adl_operationalizations.py          # 全管线
    python adl_operationalizations.py --smoke  # 仅 CHARLS + ELSA + HRS

依赖:
    - M1 抽取层 (D:\cursorproj\cursorEvaluation\cur260711\analysis\data\*.pkl)
    - CHARLS.db / HRS.db / ELSA.db / ... (SQLite)
    - smoke_charls_2pl.py 中的 2PL 算法
"""
import numpy as np
import pandas as pd
import os, sys, time, argparse, sqlite3, json, warnings
from scipy.optimize import minimize, minimize_scalar

warnings.filterwarnings('ignore')

# ======== Paths ========
HERE = os.path.dirname(os.path.abspath(__file__))
M1_DATA_DIR = r'D:\cursorproj\cursorEvaluation\cur260711\analysis\data'
OUT_DIR = os.path.join(HERE, 'analysis', 'data')
DB_DIR = r'D:\clinicdatabase\SQLitedatabase'
os.makedirs(OUT_DIR, exist_ok=True)

ADL_ITEMS = ['dressa', 'batha', 'eata', 'beda', 'toilta']
ITEM_LABELS_CN = ['穿衣', '洗澡', '进食', '上下床', '如厕']
ITEM_LABELS = ['dressing', 'bathing', 'eating', 'transfer', 'toileting']
CHARLS_4LEVEL = ['db010', 'db011', 'db012', 'db013', 'db014']

# ======== Self-implemented 2PL ========
N_ITEMS = 5
GH_NODES, GH_W = np.polynomial.hermite_e.hermegauss(21)
GH_W = GH_W / np.sqrt(2 * np.pi)

def p_2pl(theta, a, b):
    return 1.0 / (1.0 + np.exp(-a * (theta - b)))

def prev_from_latent(mu, a, b):
    th = GH_NODES + mu
    P = p_2pl(th[:, None], a[None, :], b[None, :])
    p_none = np.prod(1.0 - P, axis=1)
    return float(((1.0 - p_none) * GH_W).sum())

def loglik_bmu(params, resp, a, effect_coding=True):
    b = np.asarray(params[:N_ITEMS], float)
    mu = params[-1]
    theta = GH_NODES + mu
    P = p_2pl(theta[None, :, None], a[None, None, :], b[None, None, :])
    r = resp[:, None, :]
    logcond = (r * np.log(P + 1e-12) + (1 - r) * np.log(1 - P + 1e-12)).sum(axis=2)
    return -np.log((np.exp(logcond) * GH_W[None, :]).sum(axis=1) + 1e-300).sum()

def est_bmu(resp, a_init):
    """Estimate b (effect coding, zero-sum constraint) and mu."""
    # Rough initial values
    p = resp.mean(axis=0)
    p = np.clip(p, 0.05, 0.95)
    b_init = -np.log(p / (1 - p))  # rough item difficulty
    b_init = b_init - b_init.mean()  # zero-sum
    x0 = np.concatenate([b_init, [0.0]])

    res = minimize(loglik_bmu, x0, args=(resp, a_init, True),
                   method='Nelder-Mead',
                   options={'xatol': 1e-4, 'fatol': 1e-4, 'maxiter': 12000})
    b_hat = np.asarray(res.x[:N_ITEMS], float)
    shift = b_hat.mean()
    b_hat = b_hat - shift
    mu_hat = res.x[-1] - shift
    return b_hat, mu_hat, res.success, res.nit

def estimate_discriminations(resp, b_hat, mu_hat, a_hat):
    """Profile likelihood per item for discrimination."""
    for j in range(N_ITEMS):
        def neg_ll_a(a_j):
            a_tmp = a_hat.copy(); a_tmp[j] = float(a_j)
            b = b_hat.copy()
            theta = GH_NODES + mu_hat
            P = p_2pl(theta[None, :, None], a_tmp[None, None, :], b[None, None, :])
            r = resp[:, None, :]
            logcond = (r * np.log(P + 1e-12) + (1 - r) * np.log(1 - P + 1e-12)).sum(axis=2)
            return -np.log((np.exp(logcond) * GH_W[None, :]).sum(axis=1) + 1e-300).sum()
        res = minimize_scalar(neg_ll_a, bounds=(0.3, 5.0), method='bounded')
        a_hat[j] = res.x
    return a_hat

def iterative_2pl(resp, max_iter=10, tol=5e-4):
    """Iterative 2PL: alternate b/mu and a estimation."""
    a = np.ones(N_ITEMS)
    for it in range(max_iter):
        a_prev = a.copy()
        b, mu, ok, nit = est_bmu(resp, a)
        if not ok:
            return a, b, mu, False, it, None
        a = estimate_discriminations(resp, b, mu, a)
        delta = np.max(np.abs(a - a_prev))
        if delta < tol:
            return a, b, mu, True, it + 1, delta
    return a, b, mu, False, max_iter, delta

def test_information(theta, a, b):
    theta = np.atleast_1d(theta)
    P = p_2pl(theta[:, None], a[None, :], b[None, :])
    I = (a[None, :]**2) * P * (1 - P)
    return I.sum(axis=1), I


# ======== CHARLS 4-level extraction ========
def extract_charls_4level():
    """Extract CHARLS W4 original 4-level ADL (db010-db014) from w4_charls_CN table.
    Returns DataFrame aligned to CHARLS M1 pickle by ID."""
    con = sqlite3.connect(os.path.join(DB_DIR, 'CHARLS.db'))

    # Load M1 CHARLS pickle for ID mapping
    m1 = pd.read_pickle(os.path.join(M1_DATA_DIR, 'CHARLS.pkl'))
    # Get original IDs that correspond to M1 (CHARLS uses 'ID' column)
    # The w4 table uses 'ID' as household+person identifier
    m1_ids = set(m1['pid'].values)

    # Query w4_charls_CN for db010-db014 + ID
    cols = ['ID'] + CHARLS_4LEVEL
    q = 'SELECT ' + ','.join(f'"{c}"' for c in cols) + ' FROM "w4_charls_CN"'
    raw = pd.read_sql(q, con)
    con.close()

    df = pd.DataFrame({'pid': raw['ID']})
    for i, col in enumerate(CHARLS_4LEVEL):
        df[f'db{i:03d}_raw'] = pd.to_numeric(raw[col], errors='coerce')

    # Merge with M1 to align
    merged = m1[['pid']].merge(df, on='pid', how='left')

    # Coding: 1=no difficulty, 2=has difficulty but can do, 3=needs help, 4=cannot do
    # Missing codes: .r (refused), .d (don't know), .m (missing), .a (not asked), .s (skip)
    # These appear as NaN after pd.to_numeric

    return merged


def compute_m2_charls(charls_raw4):
    """Compute M2b (binary threshold) and M2a (severity weighted) for CHARLS.

    M2b (主分析): "需帮助或不能做 = 1" (3 or 4 -> 1)
    M2a (敏感性): 严重度加权 = mean of 5 items on 1-4 scale / 4

    Returns dict with M2b prevalence and M2a scores.
    """
    vals = charls_raw4[[f'db{i:03d}_raw' for i in range(5)]].values  # n x 5

    # M2b: strict threshold
    m2b_binary = ((vals == 3) | (vals == 4)).astype(float)
    m2b_any = (m2b_binary.sum(axis=1) >= 1).astype(float)

    # M2a: severity weighted
    valid_mask = (vals >= 1) & (vals <= 4)
    m2a_mean = np.where(valid_mask.any(axis=1),
                        np.nanmean(vals, axis=1), np.nan)
    m2a_scaled = m2a_mean / 4.0  # scale to [0, 1]

    # Count valid responses per item
    per_item = {}
    for j in range(5):
        v = vals[:, j]
        valid = (v >= 1) & (v <= 4)
        per_item[f'db{j:03d}'] = {
            'n_valid': int(valid.sum()),
            'n_1': int((v == 1).sum()),
            'n_2': int((v == 2).sum()),
            'n_3': int((v == 3).sum()),
            'n_4': int((v == 4).sum()),
            'n_missing': int((~valid).sum()),
            'm2b_pos': int(m2b_binary[valid, j].sum()),
        }

    return {
        'm2b_prev': float(m2b_any.mean()),
        'm2b_pos_count': int(m2b_any.sum()),
        'm2a_mean_score': float(np.nanmean(m2a_scaled)),
        'm2a_median_score': float(np.nanmedian(m2a_scaled)),
        'per_item': per_item,
        'm2b_binary': m2b_binary,
        'm2b_any': m2b_any,
        'm2a_score': m2a_scaled,
    }


# ======== Single-population 2PL diagnostic ========
def run_single_pop_2pl(pop_name, resp, verbose=True):
    """Run iterative 2PL on a single population and return diagnostics."""
    n = len(resp)
    a, b, mu, conv, niter, delta = iterative_2pl(resp)

    # Information diagnostics
    theta_grid = np.linspace(-3.0, 3.0, 201)
    tif, iic = test_information(theta_grid, a, b)
    mid = (theta_grid >= -1.5) & (theta_grid <= 1.5)
    sem = 1.0 / np.sqrt(np.maximum(tif, 0.01))
    rel = tif / (tif + 1)

    # Model-implied prevalence
    prev_model = prev_from_latent(mu, a, b)
    prev_obs = resp.sum(axis=1).clip(0, 1).mean()

    result = {
        'pop': pop_name,
        'n': n,
        'a': a,
        'b': b,
        'mu': mu,
        'converged': conv,
        'n_iter': niter,
        'final_delta': delta,
        'prev_obs': prev_obs,
        'prev_model': prev_model,
        'prev_diff': abs(prev_model - prev_obs),
        'TIF_at_mu': float(np.interp(mu, theta_grid, tif)),
        'TIF_mid_mean': float(tif[mid].mean()),
        'TIF_mid_min': float(tif[mid].min()),
        'rel_mid_mean': float(rel[mid].mean()),
        'rel_mid_min': float(rel[mid].min()),
        'SEM_mid_max': float(sem[mid].max()),
        'max_TIF': float(tif.max()),
        'tif': tif,
        'iic': iic,
        'theta_grid': theta_grid,
    }

    if verbose:
        status = 'CONV' if conv else f'NC(delta={delta:.4f})' if delta is not None else 'NC(maxiter)'
        print(f"  {pop_name:<15s} N={n:>6d}  mu={mu:>7.3f}  "
              f"prev_obs={prev_obs:.4f} prev_mod={prev_model:.4f}  "
              f"TIF_mid={result['TIF_mid_mean']:>5.2f}  rel={result['rel_mid_mean']:.3f}  "
              f"{status}  nit={niter}")

    return result


# ======== Main ========
def main(smoke=False):
    print("=" * 85)
    print("M2 操作化 + 13 人群 2PL 前置诊断")
    print("V9.1 Protocol - CHARLS 4-level extraction & Full-population IRT smoke")
    print("=" * 85)

    # ---- Step 1: CHARLS 4-level extraction ----
    print("\n[1] CHARLS 4-level ADL (db010-db014) extraction...")
    charls_raw = extract_charls_4level()
    m2_charls = compute_m2_charls(charls_raw)
    print(f"    CHARLS M2b (需帮助或不能=1) prevalence: {m2_charls['m2b_prev']:.4f}")
    print(f"           (n_pos = {m2_charls['m2b_pos_count']} / {len(charls_raw)})")
    print(f"    CHARLS M2a 严重度加权 mean score:    {m2_charls['m2a_mean_score']:.4f}")

    print("\n    CHARLS 4-level breakdown per item:")
    print(f"    {'Item':<20s} {'N_valid':>7s} {'1=无困难':>8s} {'2=有困难能做':>10s} "
          f"{'3=需帮助':>8s} {'4=不能做':>8s} {'缺失':>6s} {'M2b pos':>8s}")
    for j, label in enumerate(ITEM_LABELS_CN):
        d = m2_charls['per_item'][f'db{j:03d}']
        print(f"    {label:<20s} {d['n_valid']:>7d} {d['n_1']:>8d} {d['n_2']:>10d} "
              f"{d['n_3']:>8d} {d['n_4']:>8d} {d['n_missing']:>6d} {d['m2b_pos']:>8d}")

    # Compute E2 estimate
    m1_prev = pd.read_pickle(os.path.join(M1_DATA_DIR, 'CHARLS.pkl'))[ADL_ITEMS].sum(axis=1).clip(0, 1).mean()
    e2 = m1_prev - m2_charls['m2b_prev']
    print(f"\n    E2 (CHARLS): P(M1) - P(M2b) = {m1_prev:.4f} - {m2_charls['m2b_prev']:.4f} = {e2:.4f}")
    print(f"    Interpretation: switching the disability threshold from 'any difficulty'")
    print(f"    to 'needs help or cannot' reduces CHARLS ADL prevalence by {e2*100:.1f} pp.")

    # ---- Step 2: 13-population 2PL diagnostics ----
    print("\n[2] 2PL IRT diagnostics for all 13 populations...")
    populations = ['CHARLS', 'ELSA', 'HRS', 'LASI', 'MHAS',
                   'SHARE_DE', 'SHARE_CZ', 'SHARE_EE', 'SHARE_SI',
                   'SHARE_PL', 'SHARE_ES', 'SHARE_IT', 'SHARE_IL']

    if smoke:
        populations = ['CHARLS', 'ELSA', 'HRS']

    results = {}
    for pop in populations:
        pkl_path = os.path.join(M1_DATA_DIR, f'{pop}.pkl')
        if not os.path.exists(pkl_path):
            print(f"  [SKIP] {pop} - no pickle at {pkl_path}")
            continue
        df = pd.read_pickle(pkl_path)
        resp = df[ADL_ITEMS].values.astype(float)
        # Exclude rows with any NaN in ADL items (HRS has sparse missing, MHAS dressa has ~847)
        valid_rows = ~np.isnan(resp).any(axis=1)
        resp = resp[valid_rows]
        if valid_rows.sum() < len(df) * 0.85:
            print(f"  [WARN] {pop}: {len(df) - valid_rows.sum()} rows ({100*(1-valid_rows.mean()):.1f}%) dropped due to item NA")
        results[pop] = run_single_pop_2pl(pop, resp)

    # ---- Step 3: Summary table ----
    print("\n" + "=" * 85)
    print("SUMMARY: 2PL Diagnostics Across Populations")
    print("=" * 85)
    print(f"{'Pop':<15s} {'N':>6s} {'mu':>7s} {'Conv':>5s} {'nit':>4s} "
          f"{'TIF_mid':>7s} {'Rel_mid':>7s} {'SEM_max':>7s} "
          f"{'Prev_obs':>8s} {'Prev_mod':>8s} {'Diff':>6s} {'a_hat':>35s}")
    print("-" * 85)

    summary_rows = []
    for pop in populations:
        if pop not in results:
            continue
        r = results[pop]
        a_str = '/'.join([f'{v:.2f}' for v in r['a']])
        print(f"{pop:<15s} {r['n']:>6d} {r['mu']:>7.3f} {str(r['converged'])[:5]:>5s} {r['n_iter']:>4d} "
              f"{r['TIF_mid_mean']:>7.2f} {r['rel_mid_mean']:>7.3f} {r['SEM_mid_max']:>7.3f} "
              f"{r['prev_obs']:>8.4f} {r['prev_model']:>8.4f} {r['prev_diff']:>6.4f} {a_str:>35s}")
        summary_rows.append({k: r[k] for k in ['pop','n','a','b','mu','converged','n_iter',
                                                 'final_delta','prev_obs','prev_model','prev_diff',
                                                 'TIF_mid_mean','rel_mid_mean','SEM_mid_max']})

    # ---- Step 4: Cross-population comparison ----
    print("\n[3] Cross-population patterns...")
    conv_ok = sum(1 for r in results.values() if r['converged'])
    print(f"    Converged: {conv_ok}/{len(results)}")
    print(f"    Prevalence range: {min(r['prev_obs'] for r in results.values()):.3f} - "
          f"{max(r['prev_obs'] for r in results.values()):.3f}")
    print(f"    Mu range: {min(r['mu'] for r in results.values()):.3f} - "
          f"{max(r['mu'] for r in results.values()):.3f}")
    print(f"    TIF_mid range: {min(r['TIF_mid_mean'] for r in results.values()):.1f} - "
          f"{max(r['TIF_mid_mean'] for r in results.values()):.1f}")

    # SHARE internal vs cross-cohort comparison
    share_pops = [p for p in populations if p.startswith('SHARE') and p in results]
    non_share = [p for p in populations if not p.startswith('SHARE') and p in results]
    if share_pops:
        share_tif = np.mean([results[p]['TIF_mid_mean'] for p in share_pops])
        nonshare_tif = np.mean([results[p]['TIF_mid_mean'] for p in non_share])
        print(f"    SHARE mean TIF_mid: {share_tif:.1f}  vs  Non-SHARE: {nonshare_tif:.1f}")
        print(f"    (SHARE TIF is ~{share_tif/nonshare_tif*100:.0f}% of non-SHARE - "
              f"confirms small-sample caveat)")

    # ---- Save ----
    out = {
        'charls_m2': {k: v for k, v in m2_charls.items()
                      if k not in ('m2b_binary', 'm2b_any', 'm2a_score')},
        'e2_charls_pp': e2,
        'irt_results': summary_rows,
    }
    out_path = os.path.join(OUT_DIR, 'm2_operationalization_results.json')
    with open(out_path, 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=2, ensure_ascii=False, default=str)
    print(f"\nResults saved to: {out_path}")

    return results, m2_charls


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--smoke', action='store_true', help='Only CHARLS+ELSA+HRS')
    args = ap.parse_args()
    main(smoke=args.smoke)
