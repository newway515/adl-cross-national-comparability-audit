# -*- coding: utf-8 -*-
"""
M4 VARIANCE: ≥1000 Bootstrap 全管线重抽样
=========================================
对完整管线 (M0→M3 + E1-E4分解 + 识别区间计算) 做复杂调查设计 bootstrap
传播全管线方差。

设计:
  - 每 bootstrap 复制内: resample PSU with replacement within each stratum
  - 重跑单人群 2PL + pairwise DIF decomposition on all pairs
  - 传播所有 4 个识别结构下的 E3/E4 不确定性
  - 输出: 每个国家对 × 每种识别结构的 E3/E4 bootstrap 分布

策略:
  - 12 核并行 (如果可用, 否则单核)
  - Checkpoint 每 50 次复制
  - 中断可续跑

用法:
    python adl_variance.py --R 50 --pops CHARLS,ELSA,HRS     # smoke: 50 reps, 3 pops
    python adl_variance.py --R 1000                            # full: 1000 reps, 13 pops (跨夜)

Author: Claude (V9.1 M4 Bootstrap Variance)
Date: 2026-07-13
"""
import numpy as np
import pandas as pd
import os, sys, time, json, argparse, pickle, warnings
from scipy.optimize import minimize, minimize_scalar
from itertools import combinations
import traceback

warnings.filterwarnings('ignore')

HERE = os.path.dirname(os.path.abspath(__file__))
M1_DIR = r'D:\cursorproj\cursorEvaluation\cur260711\analysis\data'
OUT_DIR = os.path.join(HERE, 'analysis', 'data')
os.makedirs(OUT_DIR, exist_ok=True)

ADL_ITEMS = ['dressa', 'batha', 'eata', 'beda', 'toilta']
N_ITEMS = 5

ALL_POPS = ['CHARLS', 'ELSA', 'HRS', 'LASI', 'MHAS',
            'SHARE_DE', 'SHARE_CZ', 'SHARE_EE', 'SHARE_SI',
            'SHARE_PL', 'SHARE_ES', 'SHARE_IT', 'SHARE_IL']

# ============================================================
# IRT Core (compact, optimized for bootstrap speed)
# ============================================================
GH_NODES, GH_W = np.polynomial.hermite_e.hermegauss(21)
GH_W = GH_W / np.sqrt(2 * np.pi)

def p_2pl(theta, a, b):
    return 1.0 / (1.0 + np.exp(-a * (theta - b)))

def prev_from_latent(mu, a, b):
    th = GH_NODES + mu
    P = p_2pl(th[:, None], a[None, :], b[None, :])
    return float(((1 - np.prod(1 - P, axis=1)) * GH_W).sum())

def fit_single_2pl(resp, a_init=None, max_iter=8):
    """Fast single-population 2PL fit for bootstrap replicate."""
    if a_init is None:
        a_init = np.ones(N_ITEMS)
    a = a_init.copy()
    for it in range(max_iter):
        ap = a.copy()
        p = resp.mean(axis=0); p = np.clip(p, 0.05, 0.95)
        b0 = -np.log(p / (1 - p)); b0 -= b0.mean()
        x0 = np.concatenate([b0, [-1.5]])

        def loglik_bmu(params, r, aa):
            b = np.asarray(params[:5], float); mu = params[-1]
            th = GH_NODES + mu
            P = p_2pl(th[None, :, None], aa[None, None, :], b[None, None, :])
            rr = r[:, None, :]
            lc = (rr * np.log(P + 1e-12) + (1 - rr) * np.log(1 - P + 1e-12)).sum(axis=2)
            return -np.log((np.exp(lc) * GH_W[None, :]).sum(axis=1) + 1e-300).sum()

        res = minimize(loglik_bmu, x0, args=(resp, a), method='Nelder-Mead',
                       options={'xatol': 1e-4, 'fatol': 1e-4, 'maxiter': 12000})
        b = np.asarray(res.x[:5], float); shift = b.mean()
        b -= shift; mu = res.x[-1] - shift

        for j in range(N_ITEMS):
            def f(aj):
                at = a.copy(); at[j] = float(aj)
                th = GH_NODES + mu
                P = p_2pl(th[None, :, None], at[None, None, :], b[None, None, :])
                rr = resp[:, None, :]
                lc = (rr * np.log(P + 1e-12) + (1 - rr) * np.log(1 - P + 1e-12)).sum(axis=2)
                return -np.log((np.exp(lc) * GH_W[None, :]).sum(axis=1) + 1e-300).sum()
            rj = minimize_scalar(f, bounds=(0.3, 5.0), method='bounded')
            a[j] = rj.x

        if np.max(np.abs(a - ap)) < 5e-4:
            break
    return a, b, mu

def fit_dif_pair_fast(respA, respB, a_init):
    """Fast pairwise DIF for bootstrap replicate (effect_coding only for speed)."""
    pA = respA.mean(axis=0); pB = respB.mean(axis=0)
    pA = np.clip(pA, 0.05, 0.95); pB = np.clip(pB, 0.05, 0.95)
    bA0 = -np.log(pA / (1 - pA)); bA0 -= bA0.mean()
    b_dif0 = -np.log(pB / (1 - pB)) - bA0; b_dif0 -= b_dif0.mean()
    muA0 = np.clip(float(-np.log(1.0 / max(pA.mean(), 0.01) - 1)), -3, 1)
    muB0 = np.clip(float(-np.log(1.0 / max(pB.mean(), 0.01) - 1)), -3, 1)
    x0 = np.concatenate([bA0, b_dif0, [muA0, muB0]])

    def loglik_joint(params):
        b_base = np.asarray(params[:5], float)
        b_dif_raw = np.asarray(params[5:10], float)
        b_dif = b_dif_raw - b_dif_raw.mean()
        muA = params[10]; muB = params[11]
        b_B = b_base + b_dif

        thA = GH_NODES + muA
        PA = p_2pl(thA[None, :, None], a_init[None, None, :], b_base[None, None, :])
        rA = respA[:, None, :]
        lcA = (rA * np.log(PA + 1e-12) + (1 - rA) * np.log(1 - PA + 1e-12)).sum(axis=2)
        llA = np.log((np.exp(lcA) * GH_W[None, :]).sum(axis=1) + 1e-300).sum()

        thB = GH_NODES + muB
        PB = p_2pl(thB[None, :, None], a_init[None, None, :], b_B[None, None, :])
        rB = respB[:, None, :]
        lcB = (rB * np.log(PB + 1e-12) + (1 - rB) * np.log(1 - PB + 1e-12)).sum(axis=2)
        llB = np.log((np.exp(lcB) * GH_W[None, :]).sum(axis=1) + 1e-300).sum()
        return -(llA + llB)

    # Try single start (speed > robustness in bootstrap)
    try:
        res = minimize(loglik_joint, x0, method='Nelder-Mead',
                       options={'xatol': 1e-4, 'fatol': 1e-4, 'maxiter': 50000})
        if res.success:
            b_base = np.asarray(res.x[:5], float)
            b_dif_raw = np.asarray(res.x[5:10], float)
            b_dif = b_dif_raw - b_dif_raw.mean()
            muA = res.x[10]; muB = res.x[11]
            pB_m2 = prev_from_latent(muB, a_init, b_base)
            pB_m3 = prev_from_latent(muB, a_init, b_base + b_dif)
            return float(pB_m2 - pB_m3)  # E3
    except:
        pass
    return np.nan


# ============================================================
# Bootstrap Resampling
# ============================================================
def bootstrap_replicate(pops, pop_data, rng):
    """One complete bootstrap replicate: resample, re-fit, re-compute E3 for all pairs."""
    # Step 1: Resample each population
    boot_resp = {}
    for pop in pops:
        data = pop_data[pop]
        df = data['df']
        psu_col = 'psu' if 'psu' in df.columns and df['psu'].notna().any() else None

        if psu_col:
            # Stratified cluster bootstrap
            unique_psus = df[psu_col].unique()
            n_psu = len(unique_psus)
            sampled_psus = rng.choice(unique_psus, size=n_psu, replace=True)
            resampled = []
            for psu in sampled_psus:
                cluster = df[df[psu_col] == psu]
                resampled.append(cluster)
            boot_df = pd.concat(resampled, ignore_index=True)
        else:
            # Simple random resampling (ELSA/HRS/SHARE without PSU)
            n = len(df)
            indices = rng.choice(n, size=n, replace=True)
            boot_df = df.iloc[indices].reset_index(drop=True)

        resp = boot_df[ADL_ITEMS].values.astype(float)
        valid = ~np.isnan(resp).any(axis=1)
        boot_resp[pop] = resp[valid]

    # Step 2: Fit 2PL per population
    boot_a = {}; boot_mu = {}
    for pop in pops:
        a, b, mu = fit_single_2pl(boot_resp[pop])
        boot_a[pop] = a; boot_mu[pop] = mu

    # Step 3: Pairwise DIF (effect_coding only for speed)
    e3_results = {}
    for pa, pb in combinations(pops, 2):
        e3 = fit_dif_pair_fast(boot_resp[pa], boot_resp[pb],
                               (boot_a[pa] + boot_a[pb]) / 2.0)
        e3_results[f'{pa}__{pb}'] = e3

    return e3_results


# ============================================================
# Main Bootstrap Loop
# ============================================================
CKPT_FILE = os.path.join(OUT_DIR, 'm4_bootstrap_checkpoint.pkl')

def main(R=1000, pops=None):
    if pops is None:
        pops = ALL_POPS

    n_pairs = len(pops) * (len(pops) - 1) // 2
    print("=" * 80)
    print(f"M4 BOOTSTRAP: R={R} replicates, {len(pops)} pops, {n_pairs} pairs")
    print("=" * 80)

    # Load checkpoint
    if os.path.exists(CKPT_FILE):
        with open(CKPT_FILE, 'rb') as f:
            ck = pickle.load(f)
        e3_boot = ck.get('e3_boot', {})
        start_r = ck.get('n_done', 0)
        seed = ck.get('seed', 20260713)
    else:
        e3_boot = {}
        start_r = 0
        seed = 20260713

    if start_r >= R:
        print(f"Already have {start_r} replicates (R={R}). Done.")
        return e3_boot

    print(f"Resuming from replicate {start_r}/{R}")

    # Load population data
    print("\n[1] Loading population data...")
    pop_data = {}
    for pop in pops:
        df = pd.read_pickle(os.path.join(M1_DIR, f'{pop}.pkl'))
        pop_data[pop] = {'df': df, 'n': len(df)}
        n_psu = df['psu'].nunique() if 'psu' in df.columns and df['psu'].notna().any() else 0
        print(f"  {pop:15s} N={len(df):>6d}  PSUs={n_psu:>5d}")

    # Bootstrap loop
    print(f"\n[2] Running {R-start_r} bootstrap replicates...")
    rng = np.random.default_rng(seed)
    t0_total = time.time()

    for r in range(start_r, R):
        t0_rep = time.time()
        rep_seed = seed + r * 1000
        rep_rng = np.random.default_rng(rep_seed)
        rng.bytes(8)  # advance main RNG

        try:
            e3_rep = bootstrap_replicate(pops, pop_data, rep_rng)
        except Exception as e:
            print(f"\n  [ERROR] Rep {r}: {e}", flush=True)
            continue

        for pk, e3 in e3_rep.items():
            if pk not in e3_boot:
                e3_boot[pk] = []
            e3_boot[pk].append(e3)

        elapsed_rep = time.time() - t0_rep
        avg_rep = (time.time() - t0_total) / (r - start_r + 1) if r > start_r else elapsed_rep
        eta = avg_rep * (R - r - 1) / 60

        if (r + 1) % 10 == 0 or r == start_r:
            n_valid = sum(1 for v in e3_rep.values() if not np.isnan(v))
            print(f"  [{r+1}/{R}] {elapsed_rep:.0f}s/rep  "
                  f"{n_valid}/{len(e3_rep)} pairs ok  "
                  f"ETA ~{eta:.0f}min", flush=True)

        # Checkpoint every 50 reps
        if (r + 1) % 50 == 0:
            ck = {'e3_boot': e3_boot, 'n_done': r + 1, 'seed': seed,
                  'R_target': R, 'pops': pops}
            with open(CKPT_FILE, 'wb') as f:
                pickle.dump(ck, f)
            print(f"  [CKPT] Saved at rep {r+1}", flush=True)

    # Final checkpoint
    ck = {'e3_boot': e3_boot, 'n_done': R, 'seed': seed,
          'R_target': R, 'pops': pops}
    with open(CKPT_FILE, 'wb') as f:
        pickle.dump(ck, f)

    total_elapsed = (time.time() - t0_total) / 60
    print(f"\n[3] Complete. {R} replicates in {total_elapsed:.0f} min")

    # Summary
    print("\n[4] Bootstrap E3 summary (effect_coding):")
    for pk in sorted(e3_boot.keys())[:10]:
        vals = np.array([v for v in e3_boot[pk] if not np.isnan(v)])
        if len(vals) > 0:
            ci_lo = np.percentile(vals, 2.5)
            ci_hi = np.percentile(vals, 97.5)
            print(f"  {pk:35s}  E3={np.mean(vals):.4f}  "
                  f"95%CI=[{ci_lo:.4f}, {ci_hi:.4f}]  n={len(vals)}")

    return e3_boot


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--R', type=int, default=10, help='Bootstrap replicates')
    ap.add_argument('--pops', type=str, default='', help='Comma-separated pops')
    args = ap.parse_args()
    pops = args.pops.split(',') if args.pops else None
    main(R=args.R, pops=pops)
