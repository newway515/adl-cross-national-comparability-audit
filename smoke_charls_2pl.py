# -*- coding: utf-8 -*-
"""
CHARLS Smoke Test: V9.1 M2/M3 Pipeline First Step
==================================================
在 CHARLS 真实数据上跑 2PL IRT，验证三项前置条件：

1. 自实现 2PL MML 在真实数据上收敛
2. 5 条目 IRT 参数 (区分度 a, 难度 b) 的估计值是否合理
3. 测验信息函数 (Test Information Function) 是否在中段 (theta -1.5 to +1.5) 提供足够信息

产出：
- IRT 参数估计 (a_hat, b_hat per item)
- 条目信息曲线 (IIC) + 测验信息函数 (TIF) 图的数据
- 潜特质均值 mu_hat
- 收敛诊断 (success flag, iterations)

用法：
    python smoke_charls_2pl.py

依赖：
    - CHARLS.pkl (M1 抽取层，位于 analysis/data/)
    - C8_full_run.py 中的自实现 2PL MML 算法 (est_bmu, loglik_bmu, prev_from_latent 等)
    - numpy, scipy, pandas (myenv 已有)
    - matplotlib (仅画图用)

Author: Claude (V9.1 smoke test)
Date: 2026-07-12
"""
import numpy as np
import pandas as pd
import os, sys, time

# ---- Paths ----
HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, 'analysis', 'data')
os.makedirs(DATA_DIR, exist_ok=True)

# Import the self-implemented 2PL from C8 (if available in same dir)
# For now, inline the key functions adapted from C8_full_run.py
# to ensure standalone reproducibility

# ============================================================
# Self-implemented 2PL MML (adapted from C8_full_run.py)
# ============================================================
N_ITEMS = 5
GH_NODES, GH_W = np.polynomial.hermite_e.hermegauss(21)
GH_W = GH_W / np.sqrt(2 * np.pi)

def p_2pl(theta, a, b):
    """2PL item response probability."""
    return 1.0 / (1.0 + np.exp(-a * (theta - b)))

def prev_from_latent(mu, a, b):
    """Map latent mean mu to 'any ADL difficulty' prevalence via analytic integration."""
    th = GH_NODES + mu
    P = p_2pl(th[:, None], a[None, :], b[None, :])
    p_none = np.prod(1.0 - P, axis=1)
    return float(((1.0 - p_none) * GH_W).sum())

def loglik_bmu(params, resp, a, anchor_idx, effect_coding=False):
    """Negative log-likelihood for 2PL MML."""
    b = np.zeros(N_ITEMS)
    if effect_coding:
        b = np.asarray(params[:N_ITEMS], float)
        mu = params[-1]
    else:
        free = [j for j in range(N_ITEMS) if j not in anchor_idx]
        b[free] = params[:len(free)]
        mu = params[-1]
    theta = GH_NODES + mu
    P = p_2pl(theta[None, :, None], a[None, None, :], b[None, None, :])
    r = resp[:, None, :]
    logcond = (r * np.log(P + 1e-12) + (1 - r) * np.log(1 - P + 1e-12)).sum(axis=2)
    return -np.log((np.exp(logcond) * GH_W[None, :]).sum(axis=1) + 1e-300).sum()

def est_bmu(resp, a_init=None, anchor_idx=(4,), effect_coding=False):
    """
    Estimate item difficulties (b) and latent mean (mu) via MML.
    Returns (b_hat[5], mu_hat, success, n_iter).
    """
    from scipy.optimize import minimize

    if a_init is None:
        a_init = np.ones(N_ITEMS)  # start with all discrimination = 1

    # Initial values
    # Quick heuristic: item difficulty ~ inverse-normal of proportion correct
    p = resp.mean(axis=0)
    p = np.clip(p, 0.05, 0.95)  # avoid extreme values
    b_init = -np.sqrt(1 + a_init**2) * (np.log(p / (1 - p)))  # rough approx

    if effect_coding:
        x0 = np.concatenate([b_init, [0.0]])
        free = list(range(N_ITEMS))
    else:
        free = [j for j in range(N_ITEMS) if j not in anchor_idx]
        x0 = np.concatenate([b_init[free], [0.0]])

    res = minimize(
        loglik_bmu, x0,
        args=(resp, a_init, anchor_idx, effect_coding),
        method='Nelder-Mead',
        options={'xatol': 1e-4, 'fatol': 1e-4, 'maxiter': 12000, 'disp': False}
    )

    b_hat = np.zeros(N_ITEMS)
    if effect_coding:
        b_hat = np.asarray(res.x[:N_ITEMS], float)
        shift = b_hat.mean()
        b_hat = b_hat - shift
        mu_hat = res.x[-1] - shift
    else:
        b_hat[free] = res.x[:len(free)]
        mu_hat = res.x[-1]

    return b_hat, mu_hat, res.success, res.nit


def test_information(theta, a, b):
    """Compute test information function at points theta.
    TIF(theta) = sum_j I_j(theta) where I_j = a_j^2 * P_j * (1-P_j)"""
    theta = np.atleast_1d(theta)
    P = p_2pl(theta[:, None], a[None, :], b[None, :])
    I = (a[None, :]**2) * P * (1 - P)
    return I.sum(axis=1), I  # TIF and per-item IIC


def full_info_diagnostics(resp, a_init, b_hat, mu_hat):
    """Compute comprehensive information diagnostics."""
    # Grid for information function evaluation
    theta_grid = np.linspace(-3.0, 3.0, 201)
    tif, iic = test_information(theta_grid, a_init, b_hat)

    # Standard error of measurement: SEM(theta) = 1/sqrt(TIF(theta))
    sem = 1.0 / np.sqrt(np.maximum(tif, 0.01))

    # Reliability analogue: r(theta) = 1 - SEM^2 (when latent variance=1)
    # More precisely: r = TIF / (TIF + 1)  (for theta ~ N(0,1))
    rel = tif / (tif + 1)

    # Key summary stats
    mid_mask = (theta_grid >= -1.5) & (theta_grid <= 1.5)

    results = {
        'theta_grid': theta_grid,
        'TIF': tif,
        'IIC': iic,  # shape (n_theta, 5)
        'SEM': sem,
        'reliability': rel,
        'TIF_mid_mean': float(tif[mid_mask].mean()),
        'TIF_mid_min': float(tif[mid_mask].min()),
        'TIF_mid_max': float(tif[mid_mask].max()),
        'rel_mid_mean': float(rel[mid_mask].mean()),
        'rel_mid_min': float(rel[mid_mask].min()),
        'SEM_mid_max': float(sem[mid_mask].max()),
        'TIF_at_mean': float(np.interp(0.0, theta_grid, tif)),
        'rel_at_mean': float(np.interp(0.0, theta_grid, rel)),
        'max_info_theta': float(theta_grid[np.argmax(tif)]),
        'max_TIF': float(tif.max()),
    }
    return results


def estimate_discriminations(resp, b_hat, mu_hat, a_init):
    """
    Re-estimate discriminations by holding b and mu fixed,
    optimizing a per item via profile likelihood.

    Returns a_hat (updated).
    """
    from scipy.optimize import minimize_scalar

    a_hat = a_init.copy()

    for j in range(N_ITEMS):
        def neg_ll_a(a_j):
            a_tmp = a_hat.copy()
            a_tmp[j] = float(a_j)
            # Use effect_coding=False with anchor_idx empty and pass b as fixed
            # Reconstruct loglik directly
            b = b_hat.copy()
            theta = GH_NODES + mu_hat
            P = p_2pl(theta[None, :, None], a_tmp[None, None, :], b[None, None, :])
            r = resp[:, None, :]
            logcond = (r * np.log(P + 1e-12) + (1 - r) * np.log(1 - P + 1e-12)).sum(axis=2)
            return -np.log((np.exp(logcond) * GH_W[None, :]).sum(axis=1) + 1e-300).sum()

        res = minimize_scalar(neg_ll_a, bounds=(0.2, 3.0), method='bounded')
        a_hat[j] = res.x

    return a_hat


def iterative_2pl_estimation(resp, max_iter=5, tol=1e-3):
    """
    Iterative estimation of 2PL parameters:
    1. Start with a_init = ones(5)
    2. Estimate b and mu (holding a fixed) via effect_coding
    3. Re-estimate a (holding b and mu fixed) per item
    4. Repeat until convergence

    Returns: (a_hat, b_hat, mu_hat, converged, n_iter, history)
    """
    a_hat = np.ones(N_ITEMS)
    history = []

    for iteration in range(max_iter):
        a_prev = a_hat.copy()

        # Step 1: Estimate b, mu (holding a fixed) - use effect coding for identifiability
        b_hat, mu_hat, success, nit = est_bmu(resp, a_init=a_hat, effect_coding=True)

        if not success:
            return a_hat, b_hat, mu_hat, False, iteration, history

        # Step 2: Re-estimate a per item
        a_hat = estimate_discriminations(resp, b_hat, mu_hat, a_hat)

        # Check convergence
        delta = np.max(np.abs(a_hat - a_prev))
        prev_any = prev_from_latent(mu_hat, a_hat, b_hat)

        history.append({
            'iter': iteration,
            'a': a_hat.copy(),
            'b': b_hat.copy(),
            'mu': mu_hat,
            'prev': prev_any,
            'delta_a': delta,
        })

        if delta < tol:
            return a_hat, b_hat, mu_hat, True, iteration + 1, history

    return a_hat, b_hat, mu_hat, False, max_iter, history


# ============================================================
# CHARLS-specific analysis
# ============================================================

def load_charls():
    """Load CHARLS M1 data."""
    pkl_path = os.path.join(DATA_DIR, 'CHARLS.pkl')
    if not os.path.exists(pkl_path):
        # Try original workspace
        pkl_path = r'D:\cursorproj\cursorEvaluation\cur260711\analysis\data\CHARLS.pkl'
    df = pd.read_pickle(pkl_path)
    return df


def main():
    print("=" * 70)
    print("CHARLS Smoke Test: 2PL IRT on Real ADL Data")
    print("V9.1 Protocol - M2 Pre-flight Check")
    print("=" * 70)

    # ---- Load data ----
    print("\n[1] Loading CHARLS M1 data...")
    df = load_charls()
    ADL_ITEMS = ['dressa', 'batha', 'eata', 'beda', 'toilta']
    ITEM_LABELS = ['穿衣(dressing)', '洗澡(bathing)', '进食(eating)',
                   '上下床(transfer)', '如厕(toileting)']

    resp = df[ADL_ITEMS].values.astype(float)
    n = len(resp)
    print(f"    N = {n}")
    print(f"    Any ADL difficulty: {resp.sum(axis=1).clip(0,1).mean():.3f}")
    for j, it in enumerate(ADL_ITEMS):
        print(f"    {ITEM_LABELS[j]:20s}: pos={int(resp[:,j].sum()):5d}  "
              f"prev={resp[:,j].mean():.4f}")

    # ---- IRT Estimation ----
    print("\n[2] Iterative 2PL estimation (effect coding)...")
    t0 = time.time()
    a_hat, b_hat, mu_hat, converged, n_iter, history = iterative_2pl_estimation(resp)
    elapsed = time.time() - t0

    print(f"    Converged: {converged}")
    print(f"    Iterations: {n_iter}")
    print(f"    Time: {elapsed:.1f}s")
    print(f"    mu_hat (latent mean): {mu_hat:.4f}")
    print()
    print(f"    {'Item':<25s} {'a (disc.)':>10s} {'b (diff.)':>10s} {'Prev':>8s} {'Info @ mu':>10s}")
    print(f"    {'-'*65}")

    # Per-item diagnostics
    theta_at_mu = np.array([mu_hat])
    _, iic_at_mu = test_information(theta_at_mu, a_hat, b_hat)

    for j, label in enumerate(ITEM_LABELS):
        info = float(iic_at_mu[0, j])
        print(f"    {label:<25s} {a_hat[j]:>10.4f} {b_hat[j]:>10.4f} "
              f"{resp[:,j].mean():>8.4f} {info:>10.4f}")

    # ---- Information Diagnostics ----
    print("\n[3] Information function diagnostics...")
    info = full_info_diagnostics(resp, a_hat, b_hat, mu_hat)

    print(f"    TIF at mu={mu_hat:.2f}:        {info['TIF_at_mean']:.2f}")
    print(f"    Reliability at mu:              {info['rel_at_mean']:.3f}")
    print(f"    Max TIF (at theta={info['max_info_theta']:.2f}): {info['max_TIF']:.2f}")
    print(f"    TIF mid-range [-1.5, +1.5]:")
    print(f"        Mean:  {info['TIF_mid_mean']:.2f}")
    print(f"        Min:   {info['TIF_mid_min']:.2f}  (at worst theta in range)")
    print(f"        Max:   {info['TIF_mid_max']:.2f}")
    print(f"    Reliability mid-range [-1.5, +1.5]:")
    print(f"        Mean:  {info['rel_mid_mean']:.3f}")
    print(f"        Min:   {info['rel_mid_min']:.3f}")
    print(f"    SEM mid-range (worst):          {info['SEM_mid_max']:.3f}")

    # ---- Prevalence estimates ----
    print("\n[4] Prevalence estimates at each M-layer...")
    # M0 = M1 = Harmonized "any difficulty" prevalence (observed)
    prev_M01 = resp.sum(axis=1).clip(0, 1).mean()

    # M2b: CHARLS with higher threshold ("need help or cannot" = 3 or 4 -> 1)
    # CHARLS original 4-level: db010-db014, 1=no difficulty, 2=has difficulty but can,
    # 3=needs help, 4=cannot do
    # Harmonized maps 2/3/4 -> 1.  M2b maps 3/4 -> 1.
    # Since we only have Harmonized data in the pickle, we compute M2b by
    # re-extracting from the original CHARLS table. For now, report what we can.
    print(f"    P(M0) = P(M1) = observed 'any difficulty' prevalence: {prev_M01:.4f}")
    print(f"    Note: M2b (CHARLS 'needs help/cannot' threshold) requires")
    print(f"    original CHARLS 4-level variables (db010-db014), not available")
    print(f"    in the current pickle. This will be extracted in M2 module.")

    # M3: IRT-based prevalence from latent model
    prev_M3 = prev_from_latent(mu_hat, a_hat, b_hat)
    print(f"    P(M3) = IRT-model-implied prevalence: {prev_M3:.4f}")
    print(f"    Note: M3 is the IRT-model-implied prevalence from the 2PL.")
    print(f"    The difference P(M01) - P(M3) reflects how well the 2PL model")
    print(f"    reproduces the observed prevalence. This is NOT E3 (DIF effect)")
    print(f"    because we are looking at a single population without DIF estimation.")

    # ---- Convergence history ----
    print("\n[5] Estimation history (a parameters across iterations)...")
    print(f"    {'Iter':<6s} {'a_dress':>8s} {'a_bath':>8s} {'a_eat':>8s} "
          f"{'a_bed':>8s} {'a_toilt':>8s} {'mu':>8s} {'delta':>8s}")
    for h in history:
        a = h['a']
        print(f"    {h['iter']:<6d} {a[0]:>8.4f} {a[1]:>8.4f} {a[2]:>8.4f} "
              f"{a[3]:>8.4f} {a[4]:>8.4f} {h['mu']:>8.4f} {h['delta_a']:>8.6f}")

    # ---- Summary Assessment ----
    print("\n" + "=" * 70)
    print("SMOKE TEST SUMMARY")
    print("=" * 70)

    checks = []

    # Check 1: Convergence
    if converged:
        checks.append(("PASS", "2PL estimation converged"))
    else:
        checks.append(("WARN", f"2PL estimation did NOT converge within {n_iter} iterations"))

    # Check 2: Reasonable discrimination parameters
    a_ok = np.all((a_hat >= 0.2) & (a_hat <= 3.0))
    if a_ok:
        checks.append(("PASS", f"All a_j in [0.2, 3.0]: {a_hat.round(3)}"))
    else:
        bad = [(j, a_hat[j]) for j in range(N_ITEMS) if a_hat[j] < 0.2 or a_hat[j] > 3.0]
        checks.append(("WARN", f"Discrimination out of range: {bad}"))

    # Check 3: Mid-range information adequacy
    if info['TIF_mid_mean'] >= 4.0:
        checks.append(("PASS", f"Mean TIF in [-1.5, 1.5] = {info['TIF_mid_mean']:.1f} >= 4"))
    elif info['TIF_mid_mean'] >= 2.5:
        checks.append(("WARN", f"Mean TIF in [-1.5, 1.5] = {info['TIF_mid_mean']:.1f} (2.5-4): adequate but limited"))
    else:
        checks.append(("FAIL", f"Mean TIF in [-1.5, 1.5] = {info['TIF_mid_mean']:.1f} < 2.5: information may be insufficient for reliable DIF detection"))

    # Check 4: Mid-range reliability
    if info['rel_mid_mean'] >= 0.70:
        checks.append(("PASS", f"Mean reliability in [-1.5, 1.5] = {info['rel_mid_mean']:.3f} >= 0.70"))
    else:
        checks.append(("WARN", f"Mean reliability in [-1.5, 1.5] = {info['rel_mid_mean']:.3f} < 0.70"))

    # Check 5: Model-implied prevalence matches observed
    prev_diff = abs(prev_M3 - prev_M01)
    if prev_diff < 0.03:
        checks.append(("PASS", f"Model-implied prevalence ({prev_M3:.4f}) close to observed ({prev_M01:.4f}), diff={prev_diff:.4f}"))
    else:
        checks.append(("WARN", f"Model-implied prevalence ({prev_M3:.4f}) differs from observed ({prev_M01:.4f}) by {prev_diff:.4f}"))

    # Print all checks
    for status, msg in checks:
        symbol = {"PASS": "  [PASS]", "WARN": "  [WARN]", "FAIL": "  [FAIL]"}[status]
        print(f"{symbol}  {msg}")

    # Final verdict
    n_fail = sum(1 for s, _ in checks if s == 'FAIL')
    n_warn = sum(1 for s, _ in checks if s == 'WARN')

    print()
    if n_fail == 0 and n_warn == 0:
        print("VERDICT: All checks passed. 5-item 2PL IRT is viable for CHARLS.")
        print("Proceed to M2 operationalization and full 13-population M3 decomposition.")
    elif n_fail == 0:
        print(f"VERDICT: {n_warn} warning(s), no failures. IRT is adequate with caveats noted above.")
        print("Proceed with M2/M3 but report information limitations explicitly.")
    else:
        print(f"VERDICT: {n_fail} failure(s), {n_warn} warning(s). See issues above.")
        print("Consider whether 5-item IRT has sufficient information for DIF detection,")
        print("especially for sparse items (eating) and small-sample populations (SHARE countries).")

    # ---- Save outputs ----
    out = {
        'a_hat': a_hat,
        'b_hat': b_hat,
        'mu_hat': mu_hat,
        'converged': converged,
        'n_iter': n_iter,
        'history': history,
        'info_diag': {k: v for k, v in info.items() if k not in ('theta_grid', 'TIF', 'IIC', 'SEM', 'reliability')},
        'prev_M01': prev_M01,
        'prev_M3': prev_M3,
        'checks': [(s, m) for s, m in checks],
        'item_labels': ITEM_LABELS,
    }

    out_path = os.path.join(DATA_DIR, 'charls_smoke_results.npy')
    np.save(out_path, out, allow_pickle=True)
    print(f"\nResults saved to: {out_path}")

    return out


if __name__ == '__main__':
    main()
