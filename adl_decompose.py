# -*- coding: utf-8 -*-
"""
M3 分解模块: 年龄-性别标准化 + E1/E2/E3/E4 四结构包络

本模块执行 V9.1 方案 §7 的核心 DIF 分解管线：
  1. 年龄-性别标准化 (13 人群等权共同标准)
  2. E1/E2 点估计 (患病率量纲)
  3. E3/E4 部分识别 (四种识别结构包络)
  4. 跨人群成对差异 (78 个国家对的识别边界)

用法:
    python adl_decompose.py              # 全 13 人群 M3 分解
    python adl_decompose.py --smoke      # CHARLS+ELSA+HRS only

依赖:
    - M1抽取层 (13 人群 .pkl)
    - M2操作化结果 (CHARLS E2)
    - 自实现 2PL (est_bmu, prev_from_latent)

Author: Claude (V9.1 M3 module)
Date: 2026-07-12
"""
import numpy as np
import pandas as pd
import os, sys, time, json, argparse, warnings
from scipy.optimize import minimize, minimize_scalar
from itertools import combinations

warnings.filterwarnings('ignore')
np.set_printoptions(precision=4, suppress=True)

# ============================================================
# 0. Paths and constants
# ============================================================
HERE = os.path.dirname(os.path.abspath(__file__))
M1_DATA_DIR = r'D:\cursorproj\cursorEvaluation\cur260711\analysis\data'
OUT_DIR = os.path.join(HERE, 'analysis', 'data')
os.makedirs(OUT_DIR, exist_ok=True)

ADL_ITEMS = ['dressa', 'batha', 'eata', 'beda', 'toilta']
N_ITEMS = 5

AGE_BINS = [65, 70, 75, 80, 85, 200]
AGE_LABELS = ['65-69', '70-74', '75-79', '80-84', '85+']
SEX_LABELS = ['male', 'female']  # 1=male, 2=female

POPULATIONS = ['CHARLS', 'ELSA', 'HRS', 'LASI', 'MHAS',
               'SHARE_DE', 'SHARE_CZ', 'SHARE_EE', 'SHARE_SI',
               'SHARE_PL', 'SHARE_ES', 'SHARE_IT', 'SHARE_IL']

POP_LABELS = {
    'CHARLS': 'China (CHARLS)',
    'ELSA': 'England (ELSA)',
    'HRS': 'USA (HRS)',
    'LASI': 'India (LASI)',
    'MHAS': 'Mexico (MHAS)',
    'SHARE_DE': 'Germany (SHARE)',
    'SHARE_CZ': 'Czechia (SHARE)',
    'SHARE_EE': 'Estonia (SHARE)',
    'SHARE_SI': 'Slovenia (SHARE)',
    'SHARE_PL': 'Poland (SHARE)',
    'SHARE_ES': 'Spain (SHARE)',
    'SHARE_IT': 'Italy (SHARE)',
    'SHARE_IL': 'Israel (SHARE)',
}

# ============================================================
# 1. Self-implemented 2PL MML
# ============================================================
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
    p = resp.mean(axis=0); p = np.clip(p, 0.05, 0.95)
    b_init = -np.log(p / (1 - p))
    b_init = b_init - b_init.mean()
    x0 = np.concatenate([b_init, [0.0]])
    res = minimize(loglik_bmu, x0, args=(resp, a_init, True),
                   method='Nelder-Mead',
                   options={'xatol': 1e-4, 'fatol': 1e-4, 'maxiter': 12000})
    b = np.asarray(res.x[:N_ITEMS], float)
    shift = b.mean(); b = b - shift; mu = res.x[-1] - shift
    return b, mu, res.success, res.nit

def estimate_discriminations(resp, b_hat, mu_hat, a_hat):
    for j in range(N_ITEMS):
        def f(aj):
            at = a_hat.copy(); at[j] = float(aj)
            th = GH_NODES + mu_hat
            P = p_2pl(th[None, :, None], at[None, None, :], b_hat[None, None, :])
            r = resp[:, None, :]
            lc = (r * np.log(P + 1e-12) + (1 - r) * np.log(1 - P + 1e-12)).sum(axis=2)
            return -np.log((np.exp(lc) * GH_W[None, :]).sum(axis=1) + 1e-300).sum()
        res = minimize_scalar(f, bounds=(0.3, 5.0), method='bounded')
        a_hat[j] = res.x
    return a_hat

def iterative_2pl(resp, max_iter=10, tol=5e-4):
    a = np.ones(N_ITEMS)
    for it in range(max_iter):
        ap = a.copy()
        b, mu, ok, nit = est_bmu(resp, a)
        if not ok:
            return a, b, mu, False, it, None
        a = estimate_discriminations(resp, b, mu, a)
        delta = np.max(np.abs(a - ap))
        if delta < tol:
            return a, b, mu, True, it + 1, delta
    return a, b, mu, False, max_iter, delta


def loglik_bmu_2pop(respA, respB, a, b_free, muA, anchor_idx, effect_coding=False):
    """Joint negative log-likelihood for two populations A and B.

    Parameters:
    - respA, respB: response matrices (nA x 5, nB x 5)
    - a: discrimination parameters (5,) — assumed invariant across populations (impact-only model)
    - muA: latent mean for population A (fixed or free)
    - b_free: free difficulty parameters (subset determined by anchor_idx)
    - anchor_idx: indices of anchor items (no DIF)
    - effect_coding: if True, zero-sum constraint on DIF

    Returns: total negative log-likelihood.

    In this parameterization:
    - Population A: uses b = b_base (reference difficulties)
    - Population B: uses b = b_base + b_dif (DIF = item-level difficulty shifts)
    - muA and muB are BOTH estimated
    """
    b_base = np.zeros(N_ITEMS)

    if effect_coding:
        # b_free = [b1..b5] for population A, + [d1..d5] for DIF, + muA, muB
        # total free params = 12
        b_base = np.asarray(b_free[:N_ITEMS], float)
        b_dif = np.asarray(b_free[N_ITEMS:2*N_ITEMS], float)
        # Apply zero-sum constraint to DIF
        b_dif = b_dif - b_dif.mean()
        b_B = b_base + b_dif
        muA_val = b_free[2*N_ITEMS]
        muB_val = b_free[2*N_ITEMS + 1]
    else:
        # Non-anchor items have DIF. anchor_idx items have b_dif = 0.
        # free params: b_base for non-anchor, b_dif for non-anchor, muA, muB
        free_idx = [j for j in range(N_ITEMS) if j not in anchor_idx]
        n_free = len(free_idx)
        b_base = np.zeros(N_ITEMS)
        b_dif = np.zeros(N_ITEMS)
        for k, j in enumerate(free_idx):
            b_base[j] = b_free[k]
            b_dif[j] = b_free[k + n_free]
        b_B = b_base + b_dif
        muA_val = b_free[2*n_free]
        muB_val = b_free[2*n_free + 1]

    # Likelihood for population A
    thetaA = GH_NODES + muA_val
    PA = p_2pl(thetaA[None, :, None], a[None, None, :], b_base[None, None, :])
    rA = respA[:, None, :]
    lcA = (rA * np.log(PA + 1e-12) + (1 - rA) * np.log(1 - PA + 1e-12)).sum(axis=2)
    llA = np.log((np.exp(lcA) * GH_W[None, :]).sum(axis=1) + 1e-300).sum()

    # Likelihood for population B
    thetaB = GH_NODES + muB_val
    PB = p_2pl(thetaB[None, :, None], a[None, None, :], b_B[None, None, :])
    rB = respB[:, None, :]
    lcB = (rB * np.log(PB + 1e-12) + (1 - rB) * np.log(1 - PB + 1e-12)).sum(axis=2)
    llB = np.log((np.exp(lcB) * GH_W[None, :]).sum(axis=1) + 1e-300).sum()

    return -(llA + llB)


def est_dif_2pop(respA, respB, a_init, anchor_idx, effect_coding=False):
    """Estimate DIF between two populations.

    Returns:
        b_base_hat, b_dif_hat, muA_hat, muB_hat, success
    """
    pA = respA.mean(axis=0); pA = np.clip(pA, 0.05, 0.95)
    pB = respB.mean(axis=0); pB = np.clip(pB, 0.05, 0.95)
    bA_init = -np.log(pA / (1 - pA)); bA_init = bA_init - bA_init.mean()
    bB_init = -np.log(pB / (1 - pB))

    # crude mu estimates from prevalence
    muA_init = -np.log(1.0/max(pA.mean(), 0.01) - 1)  # rough
    muB_init = -np.log(1.0/max(pB.mean(), 0.01) - 1)

    if effect_coding:
        # bA[5], b_dif[5] (zero-sum), muA, muB = 12 params
        b_dif_init = bB_init - bA_init
        b_dif_init = b_dif_init - b_dif_init.mean()  # zero-sum
        x0 = np.concatenate([bA_init, b_dif_init, [muA_init, muB_init]])
    else:
        free_idx = [j for j in range(N_ITEMS) if j not in anchor_idx]
        n_free = len(free_idx)
        b_free_init = bA_init[free_idx]
        b_dif_init = (bB_init - bA_init)[free_idx]
        x0 = np.concatenate([b_free_init, b_dif_init, [muA_init, muB_init]])

    res = minimize(
        loglik_bmu_2pop, x0,
        args=(respA, respB, a_init, anchor_idx, effect_coding),
        method='Nelder-Mead',
        options={'xatol': 1e-4, 'fatol': 1e-4, 'maxiter': 50000}
    )

    # If Nelder-Mead fails, try again with perturbed start
    if not res.success:
        x0_perturb = x0 + np.random.default_rng(42).normal(0, 0.3, size=len(x0))
        res = minimize(
            loglik_bmu_2pop, x0_perturb,
            args=(respA, respB, a_init, anchor_idx, effect_coding),
            method='Nelder-Mead',
            options={'xatol': 1e-4, 'fatol': 1e-4, 'maxiter': 50000}
        )

    if effect_coding:
        b_base = np.asarray(res.x[:N_ITEMS], float)
        b_dif_raw = np.asarray(res.x[N_ITEMS:2*N_ITEMS], float)
        b_dif = b_dif_raw - b_dif_raw.mean()
        muA_hat = res.x[2*N_ITEMS]
        muB_hat = res.x[2*N_ITEMS + 1]
    else:
        free_idx = [j for j in range(N_ITEMS) if j not in anchor_idx]
        n_free = len(free_idx)
        b_base = np.zeros(N_ITEMS)
        b_dif = np.zeros(N_ITEMS)
        for k, j in enumerate(free_idx):
            b_base[j] = res.x[k]
            b_dif[j] = res.x[k + n_free]
        muA_hat = res.x[2*n_free]
        muB_hat = res.x[2*n_free + 1]

    return b_base, b_dif, muA_hat, muB_hat, res.success


# ============================================================
# 2. Age-Sex Standardization
# ============================================================

def standardize_prevalence(pop_df, target_n_per_stratum=None):
    """Compute age-sex standardized ADL prevalence using 13-pop equal-weight standard.

    Returns:
        pop_prev: overall standardized prevalence
        stratum_prev: dict mapping (agegrp, sex) -> prevalence
        stratum_n: dict mapping (agegrp, sex) -> count
    """
    pop_prev = 0.0
    n_strata = 0
    stratum_prev = {}
    stratum_n = {}

    for agegrp in AGE_LABELS:
        for sex_val, sex_label in [(1, 'male'), (2, 'female')]:
            mask = (pop_df['agegrp'] == agegrp) & (pop_df['sex'] == sex_val)
            subset = pop_df[mask]
            n = len(subset)
            stratum_n[(agegrp, sex_label)] = n

            if n > 0:
                any_adl = (subset[ADL_ITEMS].fillna(0).sum(axis=1) >= 1).mean()
                stratum_prev[(agegrp, sex_label)] = any_adl
                pop_prev += any_adl
                n_strata += 1
            else:
                stratum_prev[(agegrp, sex_label)] = np.nan

    return pop_prev / max(n_strata, 1), stratum_prev, stratum_n


def standardize_dif_prevalence(resp, a, b_m3, mu_m3):
    """Compute IRT-based standardized prevalence at M3 level.

    Given DIF-adjusted item parameters b_m3 and latent mean mu_m3,
    compute the model-implied prevalence via IRT integration.
    """
    return prev_from_latent(mu_m3, a, b_m3)


# ============================================================
# 3. E1/E2/E3/E4 Decomposition
# ============================================================

def decompose_single_pop(pop_name, df, a_hat, b_base, mu_hat, e2_charls=None):
    """Decompose measurement effects for a single population.

    E1 = P(M1) - P(M0). In harmonized data, M0 ~= M1, so E1 ~= 0.
    E2 = P(M2) - P(M1). Only CHARLS has E2 != 0 (16.1 pp).
    E3 = P(M3) - P(M2). DIF effect.
    E4 = E1 + E2 + E3.

    Note: P(M0) = P(M1) in harmonized data. E1 = 0 for all populations.
    P(M2) = P(M1) for non-CHARLS. For CHARLS, P(M2) is lower by E2.
    P(M3) = IRT-model-implied prevalence with b_base (pre-DIF-adjustment),
            compared to b_m3 (post-DIF-adjustment).
            In single-population mode: P(M3) = model-implied prevalence.

    Returns dict.
    """
    resp = df[ADL_ITEMS].values.astype(float)
    valid = ~np.isnan(resp).any(axis=1)
    resp = resp[valid]

    # M0 = M1 prevalence (standardized)
    p_m0_std, strata_m0, strata_n = standardize_prevalence(df)
    # crude (unstandardized) for reference
    p_m0_crude = (df[ADL_ITEMS].fillna(0).sum(axis=1) >= 1).mean()

    # E1
    e1 = 0.0

    # E2
    if pop_name == 'CHARLS' and e2_charls is not None:
        e2 = e2_charls
        p_m2 = p_m0_crude - e2  # M2 prevalence
    else:
        e2 = 0.0
        p_m2 = p_m0_crude

    # M3 prevalence (IRT-model-implied)
    p_m3_model = prev_from_latent(mu_hat, a_hat, b_base)

    # E3: in single-population, this is the gap between observed and model-implied
    # In the DIF framework: E3 = P(M2) - P(M3)
    # where P(M3) is "after removing DIF"
    # In a single population WITHOUT DIF removal, P(M3) = P(M2) (no DIF to remove)
    # So E3 is NOT computed here — it requires cross-population DIF estimation.
    # For the single-population diagnostic, we report model-implied prevalence
    # consistency as a model-fit check.

    # Actually, per V9.1: E3 = DIF displacement.
    # Single-population E3 is zero (no DIF to estimate within one population).
    e3 = 0.0  # placeholder; will be computed pairwise in M3 decomposition

    e4 = e1 + e2 + e3

    result = {
        'pop': pop_name,
        'n': len(df),
        'n_valid_irt': len(resp),
        'p_M0_std': p_m0_std,
        'p_M0_crude': p_m0_crude,
        'p_M3_model': p_m3_model,
        'E1': e1,
        'E2': e2,
        'E3': e3,
        'E4': e4,
        'a_hat': a_hat.tolist(),
        'b_base': b_base.tolist(),
        'mu_hat': float(mu_hat),
        'strata': {f'{a}_{s}': v for (a, s), v in strata_m0.items()},
        'strata_n': {f'{a}_{s}': v for (a, s), v in strata_n.items()},
    }
    return result


# ============================================================
# 4. Pairwise DIF Decomposition (Cross-Population)
# ============================================================

def pairwise_dif_decomposition(resp_dict, a_dict, pop_order, reference_pop='HRS'):
    """Compute pairwise DIF between all populations using 4 identification structures.

    For each pair (i, j):
    - Estimate 2PL DIF with 4 identification structures
    - Compute E3 displacement (DIF contribution)
    - Compute E4 total measurement displacement
    - Report the envelope (min, max) across structures

    Returns:
        pairwise_results: dict mapping (pop_i, pop_j) -> dict with E3/E4 per structure
    """
    # 4 identification structures
    structures = {
        'effect_coding': {'effect_coding': True, 'anchor_idx': None},
        'anchor_hrs': {'effect_coding': False, 'anchor_idx': (4,)},  # anchor on last item
        'sparse_anchor': {'effect_coding': False, 'anchor_idx': (3, 4)},  # anchor on last 2
        'leave_one': {'effect_coding': False, 'anchor_idx': (0, 1, 2, 3)},  # anchor on 4 (not last)
    }

    results = {}
    n_pairs = len(pop_order) * (len(pop_order) - 1) // 2
    pair_count = 0

    for i in range(len(pop_order)):
        for j in range(i + 1, len(pop_order)):
            popA = pop_order[i]
            popB = pop_order[j]
            pair_count += 1

            respA = resp_dict[popA]
            respB = resp_dict[popB]
            # Use average of single-population a_hat as initial a
            a_init = (a_dict[popA] + a_dict[popB]) / 2.0

            struct_results = {}
            for struct_name, config in structures.items():
                try:
                    b_base, b_dif, muA, muB, ok = est_dif_2pop(
                        respA, respB, a_init.copy(),
                        config['anchor_idx'],
                        config['effect_coding']
                    )

                    if not ok:
                        struct_results[struct_name] = None
                        continue

                    # Compute E3: P(M2) - P(M3,dif-adjusted)
                    # M2 prevalence in B = model prevalence with b_base (no DIF)
                    # M3 prevalence in B = model prevalence with b_base + b_dif (with DIF)
                    # But the decomposition is:
                    #   P_B^M2 = model-implied prevalence using common b (no DIF)
                    #   P_B^M3 = model-implied prevalence using b + b_dif (with DIF partial)
                    #   E3 = dif displacement = difference caused by b_dif

                    # In our framework:
                    #   b_base = difficulty after removing DIF
                    #   b_dif = DIF displacement per item
                    #   P(M2) uses b_base (common difficulty — no DIF)
                    #   P(M3) uses b_base + b_dif (population-specific — with DIF)
                    #   E3 for B relative to A: P_B^M2 - P_B^M3

                    # M2 prevalence (IRT with b_base, population-specific mu)
                    p_m2_A = prev_from_latent(muA, a_init, b_base)
                    p_m2_B = prev_from_latent(muB, a_init, b_base)

                    # M3 prevalence (IRT with b_base + b_dif)
                    b_with_dif = b_base + b_dif
                    p_m3_B = prev_from_latent(muB, a_init, b_with_dif)

                    # E3: DIF effect on prevalence
                    e3_dif = p_m2_B - p_m3_B

                    # E1 (0 in harmonized data) + E2 (only CHARLS)
                    e2_A = 0.161 if popA == 'CHARLS' else 0.0
                    e2_B = 0.161 if popB == 'CHARLS' else 0.0

                    # Total measurement displacement between A and B
                    # E4(AB) = (E1_A+E2_A+E3_A) - (E1_B+E2_B+E3_B)
                    # Since E1=0 for all:
                    e4_pair = (e2_A + e3_dif) - e2_B  # simplified

                    struct_results[struct_name] = {
                        'muA': float(muA),
                        'muB': float(muB),
                        'b_dif': b_dif.tolist(),
                        'p_m2_B': float(p_m2_B),
                        'p_m3_B': float(p_m3_B),
                        'e3_dif': float(e3_dif),
                        'e4_pair': float(e4_pair),
                    }
                except Exception as e:
                    struct_results[struct_name] = None

            results[(popA, popB)] = struct_results
            valid_structs = sum(1 for v in struct_results.values() if v is not None)
            if pair_count % 10 == 0:
                print(f"  Pair {pair_count}/{n_pairs}: {popA} vs {popB} — {valid_structs}/4 structures converged",
                      flush=True)

    return results


def compute_identification_frontier(pairwise_results, pop_order, observed_prev):
    """Compute identification frontier for all country pairs.

    For each pair:
    - Δ_obs = observed prevalence difference
    - Δ_true ∈ [Δ_obs - max(ΔE4), Δ_obs - min(ΔE4)] across structures
    - Check if 0 is in the identified set and if direction reverses.
    """
    frontier = {}
    for (popA, popB), struct_res in pairwise_results.items():
        e4_values = []
        e3_values = []
        for struct_name, res in struct_res.items():
            if res is not None:
                e4_values.append(res['e4_pair'])
                e3_values.append(res['e3_dif'])

        if len(e4_values) == 0:
            frontier[(popA, popB)] = None
            continue

        obs_diff = observed_prev[popA] - observed_prev[popB]
        e4_min = min(e4_values)
        e4_max = max(e4_values)

        # Identified set for true difference
        true_lo = obs_diff - e4_max
        true_hi = obs_diff - e4_min
        identified_width = true_hi - true_lo

        covers_zero = (true_lo <= 0 <= true_hi)
        direction_reverses = (true_lo < 0 < true_hi)  # could be either direction

        # Check if direction is preserved across all structures
        all_positive = all(v > 0 for v in e4_values) if e4_values else False
        all_negative = all(v < 0 for v in e4_values) if e4_values else False

        if true_lo > 0:
            label = 'stable'  # A > B in all structures
        elif true_hi < 0:
            label = 'stable'  # B > A in all structures
        elif covers_zero:
            label = 'unidentifiable'  # direction cannot be determined
        else:
            label = 'unidentifiable'

        frontier[(popA, popB)] = {
            'popA': popA,
            'popB': popB,
            'obs_diff': float(obs_diff),
            'E4_envelope': [float(e4_min), float(e4_max)],
            'E3_envelope': [float(min(e3_values)), float(max(e3_values))],
            'true_diff_lo': float(true_lo),
            'true_diff_hi': float(true_hi),
            'identified_width': float(identified_width),
            'covers_zero': covers_zero,
            'direction_label': label,
            'n_structures': len(e4_values),
        }

    return frontier


# ============================================================
# 5. Main pipeline
# ============================================================

def main(smoke=False):
    print("=" * 80)
    print("M3 DECOMPOSITION: 13-Population E1/E2/E3/E4 Identification Frontier")
    print("V9.1 Protocol")
    print("=" * 80)

    populations = ['CHARLS', 'ELSA', 'HRS'] if smoke else POPULATIONS

    # ---- Step 1: Load all data and run single-population 2PL ----
    print("\n[1] Loading M1 data and running single-population 2PL...")
    resp_dict = {}
    a_dict = {}
    b_dict = {}
    mu_dict = {}
    df_dict = {}
    single_results = {}

    for pop in populations:
        pkl = os.path.join(M1_DATA_DIR, f'{pop}.pkl')
        df = pd.read_pickle(pkl)
        df_dict[pop] = df

        resp = df[ADL_ITEMS].values.astype(float)
        valid = ~np.isnan(resp).any(axis=1)
        resp_valid = resp[valid]
        resp_dict[pop] = resp_valid

        a, b, mu, conv, nit, delta = iterative_2pl(resp_valid)
        a_dict[pop] = a
        b_dict[pop] = b
        mu_dict[pop] = mu

        # Single-population decomposition
        res = decompose_single_pop(pop, df, a, b, mu, e2_charls=0.161)
        single_results[pop] = res
        print(f"  {pop:15s} N={res['n']:>6d}  mu={mu:>7.3f}  "
              f"P(M0)={res['p_M0_crude']:.4f}  P(M3)={res['p_M3_model']:.4f}  "
              f"E2={res['E2']:.4f}  conv={conv} nit={nit}",
              flush=True)

    # ---- Step 2: Pairwise DIF decomposition ----
    print(f"\n[2] Pairwise DIF decomposition ({len(populations)} pops → "
          f"{len(populations)*(len(populations)-1)//2} pairs)...")

    pairwise = pairwise_dif_decomposition(resp_dict, a_dict, populations)

    # ---- Step 3: Identification frontier ----
    print("\n[3] Computing identification frontier...")
    observed_prev = {pop: single_results[pop]['p_M0_crude'] for pop in populations}
    frontier = compute_identification_frontier(pairwise, populations, observed_prev)

    # Summary counts
    n_stable = sum(1 for v in frontier.values() if v and v['direction_label'] == 'stable')
    n_unident = sum(1 for v in frontier.values() if v and v['direction_label'] == 'unidentifiable')
    n_total = sum(1 for v in frontier.values() if v is not None)
    print(f"    Frontier: {n_stable}/{n_total} stable, {n_unident}/{n_total} unidentifiable")

    if n_total > 0:
        widths = [v['identified_width'] for v in frontier.values() if v]
        print(f"    Width range: {min(widths):.3f} – {max(widths):.3f}")
        print(f"    Median width: {np.median(widths):.3f}")
    else:
        widths = []

    # ---- Step 4: Output ----
    print(f"\n[4] Generating output tables...")

    # Guard against zero pairs
    if n_total == 0:
        print("    WARNING: No valid pairwise DIF results. Cannot compute frontier.")
        out = {
            'single_population': rows,
            'frontier_summary': {'n_pairs': 0, 'n_stable': 0, 'n_unidentifiable': 0,
                                'width_range': [], 'width_median': None},
            'top_uncertain_pairs': [],
        }
        out_path = os.path.join(OUT_DIR, 'm3_decomposition_results.json')
        with open(out_path, 'w', encoding='utf-8') as f:
            json.dump(out, f, indent=2, ensure_ascii=False, default=str)
        print(f"Results saved to: {out_path}")
        print("\nM3 point-estimate phase ABORTED — pairwise DIF estimation needs debugging.")
        print("Suggestion: run with fewer populations or increase maxiter.")
        return out

    # Table of single-population results
    rows = []
    for pop in populations:
        r = single_results[pop]
        rows.append({
            'pop': pop,
            'n': r['n'],
            'p_M0_crude': round(r['p_M0_crude'], 4),
            'p_M0_std': round(r['p_M0_std'], 4),
            'p_M3_model': round(r['p_M3_model'], 4),
            'mu_hat': round(r['mu_hat'], 4),
            'E1': round(r['E1'], 4),
            'E2': round(r['E2'], 4),
            'a_hat': [round(x, 3) for x in r['a_hat']],
            'b_base': [round(x, 4) for x in r['b_base']],
        })

    # Top-N problematic pairs (largest E4 envelope)
    frontier_sorted = sorted(
        [(k, v) for k, v in frontier.items() if v],
        key=lambda x: x[1]['identified_width'],
        reverse=True
    )

    top_pairs = []
    for (popA, popB), v in frontier_sorted[:10]:
        top_pairs.append({
            'pair': f'{popA} vs {popB}',
            'obs_diff_pp': round(v['obs_diff'] * 100, 1),
            'E4_envelope_pp': [round(x * 100, 1) for x in v['E4_envelope']],
            'true_diff_range_pp': [round(v['true_diff_lo'] * 100, 1),
                                   round(v['true_diff_hi'] * 100, 1)],
            'direction': v['direction_label'],
            'width_pp': round(v['identified_width'] * 100, 1),
        })

    out = {
        'single_population': rows,
        'frontier_summary': {
            'n_pairs': n_total,
            'n_stable': n_stable,
            'n_unidentifiable': n_unident,
            'width_range': [float(min(widths)), float(max(widths))] if widths else [],
            'width_median': float(np.median(widths)) if widths else None,
        },
        'top_uncertain_pairs': top_pairs,
        # Full pairwise results (large, for downstream use)
        '_pairwise': {f'{k[0]}__{k[1]}': v for k, v in pairwise.items()},
        '_frontier': {f'{k[0]}__{k[1]}': v for k, v in frontier.items()},
    }

    out_path = os.path.join(OUT_DIR, 'm3_decomposition_results.json')
    with open(out_path, 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=2, ensure_ascii=False, default=str)
    print(f"\nResults saved to: {out_path}")

    # ---- Step 5: Summary printout ----
    print("\n" + "=" * 80)
    print("M3 DECOMPOSITION SUMMARY")
    print("=" * 80)
    print(f"Populations: {len(populations)}")
    print(f"Country pairs: {n_total}")
    print(f"Identification frontier stable: {n_stable}/{n_total} ({100*n_stable/n_total:.0f}%)")
    print(f"Identification frontier unidentifiable: {n_unident}/{n_total} ({100*n_unident/n_total:.0f}%)")
    print()

    if top_pairs:
        print("Top uncertain pairs (widest identification frontier):")
        print(f"{'Pair':<30s} {'Obs diff':>8s} {'E4 range':>14s} {'True range':>16s} {'Width':>8s} {'Label':>16s}")
        print("-" * 95)
        for p in top_pairs[:8]:
            print(f"{p['pair']:<30s} {p['obs_diff_pp']:>7.1f}pp "
                  f"[{p['E4_envelope_pp'][0]:>5.1f}, {p['E4_envelope_pp'][1]:>5.1f}]pp "
                  f"[{p['true_diff_range_pp'][0]:>5.1f}, {p['true_diff_range_pp'][1]:>5.1f}]pp "
                  f"{p['width_pp']:>7.1f}pp {p['direction']:>16s}")

    print("\n" + "=" * 80)
    print("M3 point-estimate phase complete.")
    print("Next: M4 variance (≥1000 bootstrap) + nonparametric DIF + LD diagnostics")
    print("=" * 80)

    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--smoke', action='store_true', help='CHARLS+ELSA+HRS only')
    args = ap.parse_args()
    main(smoke=args.smoke)
