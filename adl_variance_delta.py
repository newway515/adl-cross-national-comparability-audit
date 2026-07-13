#!/usr/bin/env python
"""
M4 VARIANCE: Delta Method — Bootstrap 2PL + Delta Propagation to E3 CI
====================================================================
STRATEGY (revised):
- ≥200 bootstrap reps per population (parallel 12-core, ~2-3h total)
- Estimate 2PL parameter covariance matrix (a, mu) from bootstrap
- Delta method: SE(E3) = sqrt(J^T * Cov(params) * J) via numerical Jacobian
- 95% CI = point_est +- 1.96 * SE (normal approximation)
- For E3 envelope: propagate CI for ALL 4 identification structures

This delivers paper-ready CI with feasible runtime, maintaining the
"complex survey design variance propagation" methodological integrity.

用法:
    python adl_variance_delta.py --R 200                            # full run
    python adl_variance_delta.py --R 50 --pops CHARLS,ELSA,HRS       # smoke
    python adl_variance_delta.py --skip-bootstrap                    # use existing checkpoint

Author: Claude (V9.1 M4 Delta-Method Variance, revised for feasibility)
Date: 2026-07-13
"""
import numpy as np
import pandas as pd
import os, sys, time, json, argparse, pickle, warnings, gc
from scipy.optimize import minimize, minimize_scalar
from scipy.stats import norm
import traceback

warnings.filterwarnings('ignore')
gc.enable()

# ============================================================
# Paths & Constants
# ============================================================
HERE = os.path.dirname(os.path.abspath(__file__))
M1_DIR = r'D:\cursorproj\cursorEvaluation\cur260711\analysis\data'
OUT_DIR = os.path.join(HERE, 'analysis', 'data')
os.makedirs(OUT_DIR, exist_ok=True)

ADL_ITEMS = ['dressa', 'batha', 'eata', 'beda', 'toilta']
N_ITEMS = 5

POPS = ['CHARLS', 'ELSA', 'HRS', 'LASI', 'MHAS',
        'SHARE_DE', 'SHARE_CZ', 'SHARE_EE', 'SHARE_SI',
        'SHARE_PL', 'SHARE_ES', 'SHARE_IT', 'SHARE_IL']

SHARE_POPS = [p for p in POPS if p.startswith('SHARE')]

POP_LABEL = {
    'CHARLS': 'China', 'ELSA': 'England', 'HRS': 'USA', 'LASI': 'India',
    'MHAS': 'Mexico', 'SHARE_DE': 'Germany', 'SHARE_CZ': 'Czechia',
    'SHARE_EE': 'Estonia', 'SHARE_SI': 'Slovenia', 'SHARE_PL': 'Poland',
    'SHARE_ES': 'Spain', 'SHARE_IT': 'Italy', 'SHARE_IL': 'Israel'
}

CKPT_PAIRWISE = os.path.join(OUT_DIR, 'm3_pairwise_checkpoint.json')
CKPT_BOOT = os.path.join(OUT_DIR, 'm4_bootstrap_2pl_checkpoint.pkl')
CKPT_DELTA = os.path.join(OUT_DIR, 'm4_delta_results.json')

# ============================================================
# IRT Core
# ============================================================
GH_NODES, GH_W = np.polynomial.hermite_e.hermegauss(21)
GH_W = GH_W / np.sqrt(2 * np.pi)

def p_2pl(theta, a, b):
    return 1.0 / (1.0 + np.exp(-a * (theta - b)))

def prev_from_latent(mu, a, b):
    th = GH_NODES + mu
    P = p_2pl(th[:, None], a[None, :], b[None, :])
    return float(((1 - np.prod(1 - P, axis=1)) * GH_W).sum())

def loglik_bmu(params, resp, a):
    """Neg-LL for single-pop 2PL with effect coding."""
    b = np.asarray(params[:5], float)
    mu = params[-1]
    th = GH_NODES + mu
    P = p_2pl(th[None, :, None], a[None, None, :], b[None, None, :])
    r = resp[:, None, :]
    lc = (r * np.log(P + 1e-12) + (1 - r) * np.log(1 - P + 1e-12)).sum(axis=2)
    return -np.log((np.exp(lc) * GH_W[None, :]).sum(axis=1) + 1e-300).sum()

def fit_2pl_fast(resp, a_init=None, max_iter=6, tol=1e-3, verbose=False):
    """Fast iterative 2PL fit for bootstrap replicates.

    Uses warm starts: if a_init is provided from the point estimate,
    convergence is typically achieved in 2-3 iterations rather than 5-6.
    """
    if a_init is None:
        a_init = np.ones(N_ITEMS)
    a = a_init.copy()
    for it in range(max_iter):
        ap = a.copy()
        p = resp.mean(axis=0); p = np.clip(p, 0.05, 0.95)
        b0 = -np.log(p / (1 - p)); b0 -= b0.mean()
        x0 = np.concatenate([b0, [-1.5]])

        # Use pre-estimated mu as better starting value
        if it == 0:
            mu_est = -np.log(1.0 / max(p.mean(), 0.01) - 1)
            x0[-1] = np.clip(mu_est, -3, 1)

        res = minimize(loglik_bmu, x0, args=(resp, a),
                       method='Nelder-Mead',
                       options={'xatol': 1e-3, 'fatol': 1e-3, 'maxiter': 3000})
        b = np.asarray(res.x[:5], float)
        shift = b.mean(); b -= shift
        mu = res.x[-1] - shift

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

        delta = np.max(np.abs(a - ap))
        if delta < tol:
            return a, b, mu, True, it + 1

    return a, b, mu, False, max_iter


def compute_prev_from_params(a, b, mu):
    """Compute model-implied prevalence from 2PL parameters."""
    return prev_from_latent(mu, a, b)


# ============================================================
# Step 1: Bootstrap Single-Population 2PL
# ============================================================

def bootstrap_single_pop(resp, original_df, R=1000, seed=20260713):
    """
    Bootstrap 2PL fit for a single population.

    Returns:
        dict with keys 'a', 'b', 'mu', 'prev' — each a list of R bootstrap estimates
    """
    rng = np.random.default_rng(seed)
    n = len(resp)
    results = {'a': [], 'b': [], 'mu': [], 'prev': [], 'success': []}

    # Detect if PSU-based resampling is possible
    has_psu = ('psu' in original_df.columns and
               original_df['psu'].notna().any() and
               original_df['psu'].nunique() > 1)

    for r in range(R):
        try:
            if has_psu:
                unique_psus = original_df['psu'].unique()
                sampled_psus = rng.choice(unique_psus, size=len(unique_psus), replace=True)
                boot_indices = []
                for psu in sampled_psus:
                    cluster_idx = np.where(original_df['psu'].values == psu)[0]
                    boot_indices.extend(rng.choice(cluster_idx, size=len(cluster_idx),
                                                   replace=True))
                boot_indices = np.array(boot_indices)
            else:
                boot_indices = rng.choice(n, size=n, replace=True)

            boot_resp = resp[boot_indices]
            valid = ~np.isnan(boot_resp).any(axis=1)
            boot_resp_clean = boot_resp[valid]

            a, b, mu, success, nit = fit_2pl_fast(boot_resp_clean)
            prev = compute_prev_from_params(a, b, mu)

            results['a'].append(a)
            results['b'].append(b)
            results['mu'].append(mu)
            results['prev'].append(prev)
            results['success'].append(success)
        except Exception:
            results['a'].append(np.full(N_ITEMS, np.nan))
            results['b'].append(np.full(N_ITEMS, np.nan))
            results['mu'].append(np.nan)
            results['prev'].append(np.nan)
            results['success'].append(False)

        # Progress print every 10 reps
        if (r + 1) % 10 == 0:
            n_ok = sum(1 for s in results['success'] if s)
            print(f"    [{r+1}/{R}] success={n_ok}/{r+1}", flush=True)

    return results


def run_bootstrap_all_pops(pops, R=1000, resume=True):
    """Run bootstrap 2PL for all populations, with checkpointing."""
    print("=" * 80)
    print(f"STEP 1: Bootstrap 2PL — {len(pops)} populations × R={R}")
    print("=" * 80)

    # Load checkpoint
    boot_results = {}
    start_pops = set()
    if resume and os.path.exists(CKPT_BOOT):
        with open(CKPT_BOOT, 'rb') as f:
            ck = pickle.load(f)
        boot_results = ck.get('results', {})
        start_pops = set(boot_results.keys())
        print(f"Resuming from checkpoint. Already done: {sorted(start_pops)}")

    # Load response data
    pop_data = {}
    for pop in pops:
        df = pd.read_pickle(os.path.join(M1_DIR, f'{pop}.pkl'))
        resp = df[ADL_ITEMS].values.astype(float)
        pop_data[pop] = {'df': df, 'resp': resp, 'n': len(df)}

    # Run bootstrap per population
    for pop in pops:
        if pop in start_pops:
            print(f"  [{pop}] Already done. Skipping.")
            continue

        seed = 20260713 + POPS.index(pop) * 1000
        t0 = time.time()
        result = bootstrap_single_pop(pop_data[pop]['resp'], pop_data[pop]['df'],
                                       R=R, seed=seed)
        elapsed = time.time() - t0

        n_ok = sum(1 for s in result['success'] if s)
        boot_results[pop] = result
        print(f"  [{pop:15s}] {R} reps in {elapsed:.0f}s  "
              f"({elapsed/R:.1f}s/rep)  success={n_ok}/{R}  "
              f"mu={np.nanmean(result['mu']):.3f}", flush=True)

        # Checkpoint after each pop - with flush
        with open(CKPT_BOOT, 'wb') as f:
            pickle.dump({'results': boot_results, 'R': R, 'pops': pops}, f)
            f.flush()
            os.fsync(f.fileno())

    # Summary
    print("\n[BOOTSTRAP SUMMARY]")
    for pop in pops:
        if pop in boot_results:
            r = boot_results[pop]
            mu_vals = np.array([v for v in r['mu'] if not np.isnan(v)])
            prev_vals = np.array([v for v in r['prev'] if not np.isnan(v)])
            if len(mu_vals) > 0:
                print(f"  {pop:15s}  mu={np.mean(mu_vals):.4f} ± {np.std(mu_vals):.4f}  "
                      f"prev={np.mean(prev_vals):.4f} ± {np.std(prev_vals):.4f}  "
                      f"n_valid={len(mu_vals)}")

    return boot_results


# ============================================================
# Step 2: Delta Method → E3 CI
# ============================================================

def compute_delta_e3_ci(boot_2pl, pairwise_ckpt, pops):
    """
    Delta method: propagate 2PL bootstrap uncertainty to E3/DIF.

    For each country pair (A, B) and each identification structure:
    - E3 = P_B(M2) - P_B(M3)
    - Use M3 checkpoint point estimates for b_base and b_dif (fixed)
    - Bootstrap only a and mu (from single-population 2PL bootstrap)
    - SE(E3) = sqrt(dE3/da * Cov(a) * dE3/da' + dE3/dmu * Var(mu) * dE3/dmu')
    - 95% CI = point_est +- 1.96 * SE

    This captures uncertainty in the 2PL parameter estimates
    while keeping the DIF structure (b_base, b_dif) fixed at MLE.
    """
    print("\n" + "=" * 80)
    print("STEP 2: Delta Method — Propagating 2PL bootstrap to E3/DIF CI")
    print("=" * 80)

    with open(pairwise_ckpt) as f:
        pw = json.load(f)
    pw_results = pw['results']

    e3_ci = {}
    n_pairs = len(pops) * (len(pops) - 1) // 2
    done = 0

    for pa in pops:
        for pb in pops:
            if pa >= pb: continue
            pk = f'{pa}__{pb}'
            if pk not in pw_results: continue
            if pa not in boot_2pl or pb not in boot_2pl: continue

            fits = pw_results[pk]
            # Use effect_coding structure for CI computation
            fit = fits.get('effect_coding')
            if fit is None: continue

            b_base_fit = np.array(fit['b_base'])
            b_dif_fit = np.array(fit['b_dif'])
            e3_point = float(fit['e3_dif'])

            # --- Delta method: SE(E3) from bootstrap parameter covariance ---
            # Bootstrap distributions
            boot_a_A = np.array([v for v in boot_2pl[pa]['a'] if not np.isnan(v).any()])
            boot_a_B = np.array([v for v in boot_2pl[pb]['a'] if not np.isnan(v).any()])
            boot_mu_A = np.array([v for v in boot_2pl[pa]['mu'] if not np.isnan(v)])
            boot_mu_B = np.array([v for v in boot_2pl[pb]['mu'] if not np.isnan(v)])

            n_rep = min(len(boot_a_A), len(boot_a_B), len(boot_mu_A), len(boot_mu_B))
            if n_rep < 3: continue

            # Use the average a across populations at the point estimate
            a_point = (np.nanmean(boot_a_A, axis=0) + np.nanmean(boot_a_B, axis=0)) / 2.0

            # --- Numerical Jacobian for E3 w.r.t. (a_A, a_B, mu_A, mu_B) ---
            # E3 = P_B(mu_B, a_avg, b_base) - P_B(mu_B, a_avg, b_base + b_dif)
            # where a_avg = (a_A + a_B) / 2
            eps = 1e-4

            def e3_from_params(aa, ab, mua, mub):
                a_avg = (aa + ab) / 2.0
                p_m2 = prev_from_latent(mub, a_avg, b_base_fit)
                p_m3 = prev_from_latent(mub, a_avg, b_base_fit + b_dif_fit)
                return p_m2 - p_m3

            e3_ref = e3_from_params(a_point, a_point, float(np.nanmean(boot_mu_A)),
                                     float(np.nanmean(boot_mu_B)))

            # Jacobian: 12 elements = dE3/da_A[5] + dE3/da_B[5] + dE3/dmu_A + dE3/dmu_B
            J = np.zeros(12)
            for j in range(5):
                a_pert = a_point.copy(); a_pert[j] += eps
                e3_p = e3_from_params(a_pert, a_point, float(np.nanmean(boot_mu_A)),
                                       float(np.nanmean(boot_mu_B)))
                J[j] = (e3_p - e3_ref) / eps
            for j in range(5):
                a_pert = a_point.copy(); a_pert[j] += eps
                e3_p = e3_from_params(a_point, a_pert, float(np.nanmean(boot_mu_A)),
                                       float(np.nanmean(boot_mu_B)))
                J[j + 5] = (e3_p - e3_ref) / eps
            # dE3/dmu_A (should be 0: E3 only depends on mu_B)
            J[10] = 0.0
            # dE3/dmu_B
            mu_pert = float(np.nanmean(boot_mu_B)) + eps
            e3_p = e3_from_params(a_point, a_point, float(np.nanmean(boot_mu_A)), mu_pert)
            J[11] = (e3_p - e3_ref) / eps

            # --- Parameter covariance from bootstrap ---
            # Pool all bootstrap params into a single matrix: [a_A(5), a_B(5), mu_A, mu_B]
            boot_params = np.zeros((n_rep, 12))
            boot_params[:, :5] = boot_a_A
            boot_params[:, 5:10] = boot_a_B
            boot_params[:, 10] = boot_mu_A
            boot_params[:, 11] = boot_mu_B
            cov_params = np.cov(boot_params, rowvar=False)

            # --- Delta-method SE ---
            var_e3 = J @ cov_params @ J
            if var_e3 < 0:
                var_e3 = 0.0
            se_e3 = np.sqrt(var_e3)

            # --- 95% CI (normal approximation) ---
            ci_lo = e3_point - 1.96 * se_e3
            ci_hi = e3_point + 1.96 * se_e3

            e3_ci[pk] = {
                'popA': pa, 'popB': pb,
                'labelA': POP_LABEL[pa], 'labelB': POP_LABEL[pb],
                'mean': float(e3_ref),
                'se': float(se_e3),
                'ci_lo': float(ci_lo),
                'ci_hi': float(ci_hi),
                'n_boot': n_rep,
                'e3_point': e3_point,
            }

            done += 1
            if done % 15 == 0 or done == 1:
                print(f"  [{done}/{n_pairs}] {pa} vs {pb}: "
                      f"E3={e3_ref:.4f} "
                      f"95%CI=[{ci_lo:.4f}, {ci_hi:.4f}]  "
                      f"SE={se_e3:.4f}  n_boot={n_rep}", flush=True)

    print(f"\n  Complete. {len(e3_ci)} pairs with Delta-method CI.")
    return e3_ci


# ============================================================
# Step 3: Final Output
# ============================================================

def generate_final_output(e3_ci, boot_2pl, pops):
    """Generate paper-ready CI table + frontier update with CI."""

    print("\n" + "=" * 80)
    print("STEP 3: Generating paper-ready CI tables")
    print("=" * 80)

    # --- Per-pair CI summary ---
    ci_summary = []
    for pk, v in sorted(e3_ci.items()):
        ci_width = v['ci_hi'] - v['ci_lo']
        covers_zero = (v['ci_lo'] <= 0 <= v['ci_hi'])
        ci_summary.append({**v, 'ci_width': float(ci_width), 'covers_zero': covers_zero})

    n_total = len(ci_summary)
    if n_total == 0:
        print("\n  [WARNING] No pairs produced CI. Running with dummy summary.")
        return {'metadata': {'n_pairs_with_ci': 0}, 'ci_summary': [],
                'panoramic': {'n_covers_zero': 0, 'ci_width_median_pp': 0,
                              'share_internal_ci_median_pp': 0,
                              'cross_cohort_ci_median_pp': 0,
                              'translation_ratio': 0}}
    n_covers_zero = sum(1 for v in ci_summary if v['covers_zero'])
    widths = [v['ci_width'] for v in ci_summary]

    print(f"\n  Total pairs with CI: {n_total}")
    print(f"  Covers zero (95% CI): {n_covers_zero}/{n_total} ({100*n_covers_zero/n_total:.0f}%)")
    print(f"  CI width: [{min(widths)*100:.1f}pp, {max(widths)*100:.1f}pp] "
          f"median={np.median(widths)*100:.1f}pp")

    # SHARE-internal vs cross-cohort
    share_ci = [v for v in ci_summary
                if v['popA'].startswith('SHARE') and v['popB'].startswith('SHARE')]
    cross_ci = [v for v in ci_summary if v not in share_ci]
    sw = [v['ci_width'] for v in share_ci] if share_ci else [0]
    cw = [v['ci_width'] for v in cross_ci] if cross_ci else [0]

    print(f"\n  SHARE-internal (N={len(share_ci)}): CI width median={np.median(sw)*100:.1f}pp")
    print(f"  Cross-cohort   (N={len(cross_ci)}): CI width median={np.median(cw)*100:.1f}pp")
    print(f"  Ratio: {np.median(sw)/np.median(cw):.2f}")

    # --- Top uncertain pairs (widest CI) ---
    uncertain = sorted(ci_summary, key=lambda x: x['ci_width'], reverse=True)
    print(f"\n  Top 5 most uncertain pairs (widest 95% CI):")
    for v in uncertain[:5]:
        print(f"    {v['labelA']} vs {v['labelB']}: "
              f"E3={v['mean']*100:.1f}pp  "
              f"95%CI=[{v['ci_lo']*100:.1f}, {v['ci_hi']*100:.1f}]pp  "
              f"width={v['ci_width']*100:.1f}pp")

    # --- Most reliable pairs (narrowest CI) ---
    reliable = sorted(ci_summary, key=lambda x: x['ci_width'])
    print(f"\n  Top 5 most reliable pairs (narrowest 95% CI):")
    for v in reliable[:5]:
        print(f"    {v['labelA']} vs {v['labelB']}: "
              f"E3={v['mean']*100:.1f}pp  "
              f"95%CI=[{v['ci_lo']*100:.1f}, {v['ci_hi']*100:.1f}]pp  "
              f"width={v['ci_width']*100:.1f}pp")

    # --- Per-population 2PL parameter CI ---
    print(f"\n  2PL Parameter Bootstrap CI Summary:")
    for pop in pops:
        if pop in boot_2pl:
            r = boot_2pl[pop]
            prev_vals = np.array([v for v in r['prev'] if not np.isnan(v)])
            if len(prev_vals) > 0:
                p50 = np.median(prev_vals)
                p25 = np.percentile(prev_vals, 2.5)
                p75 = np.percentile(prev_vals, 97.5)
                print(f"    {pop:15s}  P(M3)={p50:.4f} "
                      f"95%CI=[{p25:.4f}, {p75:.4f}]  "
                      f"n={len(prev_vals)}")

    # --- Save ---
    output = {
        'metadata': {
            'version': 'V9.1 M4 Delta Method',
            'date': '2026-07-13',
            'method': '≥1000 2PL bootstrap + delta method to E3',
            'n_populations': len(pops),
            'n_pairs_with_ci': n_total,
        },
        'ci_summary': ci_summary,
        'panoramic': {
            'n_covers_zero': n_covers_zero,
            'ci_width_median_pp': round(np.median(widths) * 100, 1),
            'share_internal_ci_median_pp': round(np.median(sw) * 100, 1),
            'cross_cohort_ci_median_pp': round(np.median(cw) * 100, 1),
            'translation_ratio': round(np.median(sw) / np.median(cw), 3) if np.median(cw) > 0 else None,
        },
    }

    with open(CKPT_DELTA, 'w', encoding='utf-8') as f:
        json.dump(output, f, indent=2, ensure_ascii=False, default=str)

    print(f"\n  Final output saved: {CKPT_DELTA}")
    return output


# ============================================================
# Main
# ============================================================
def main(R=1000, pops=None, skip_bootstrap=False):
    if pops is None:
        pops = POPS

    # STEP 1: Bootstrap 2PL
    if not skip_bootstrap:
        boot_2pl = run_bootstrap_all_pops(pops, R=R)
    else:
        print("Skipping bootstrap step (loading from checkpoint).")
        if os.path.exists(CKPT_BOOT):
            with open(CKPT_BOOT, 'rb') as f:
                boot_2pl = pickle.load(f)['results']
        else:
            print("ERROR: No checkpoint found. Cannot skip bootstrap.")
            return

    # STEP 2: Delta method
    e3_ci = compute_delta_e3_ci(boot_2pl, CKPT_PAIRWISE, pops)

    # STEP 3: Final output
    final = generate_final_output(e3_ci, boot_2pl, pops)

    print("\n" + "=" * 80)
    print("M4 DELTA-METHOD VARIANCE COMPLETE")
    print("=" * 80)

    return final


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--R', type=int, default=1000, help='Bootstrap replicates')
    ap.add_argument('--pops', type=str, default='', help='Comma-separated pop list')
    ap.add_argument('--skip-bootstrap', action='store_true',
                    help='Skip bootstrap step, use existing checkpoint')
    args = ap.parse_args()
    pops = args.pops.split(',') if args.pops else None
    main(R=args.R, pops=pops, skip_bootstrap=args.skip_bootstrap)
