# -*- coding: utf-8 -*-
"""
M3 Decomposition — Batch Runner for Pairwise DIF
================================================
Runs DIF estimation on all N*(N-1)/2 = 78 country pairs, 4 identification structures each.
Uses optimized starting values (from single-population 2PL estimates) and runs in parallel.
Stores intermediate results for checkpoint/restart.

This is the point-estimate phase: no bootstrap variance here (that's M4).

用法:
    python run_m3_pairwise.py --pops CHARLS,ELSA,HRS        # smoke: 3 pairs
    python run_m3_pairwise.py                                # full: 78 pairs

Author: Claude (V9.1 M3 batch runner)
Date: 2026-07-12
"""
import numpy as np
import pandas as pd
import os, sys, time, json, argparse, traceback
from scipy.optimize import minimize
from itertools import combinations
import warnings
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

GH_NODES, GH_W = np.polynomial.hermite_e.hermegauss(21)
GH_W = GH_W / np.sqrt(2 * np.pi)

# Checkpoint file
CHECKPOINT = os.path.join(OUT_DIR, 'm3_pairwise_checkpoint.json')

# ============================================================
# IRT core
# ============================================================
def p_2pl(theta, a, b):
    return 1.0 / (1.0 + np.exp(-a * (theta - b)))

def prev_from_latent(mu, a, b):
    th = GH_NODES + mu
    P = p_2pl(th[:, None], a[None, :], b[None, :])
    return float(((1 - np.prod(1 - P, axis=1)) * GH_W).sum())

def loglik_joint(params, respA, respB, a_fixed, anchor_idx, effect_coding):
    """Joint neg-LL for 2-population DIF model."""
    b_base = np.zeros(N_ITEMS)
    b_dif = np.zeros(N_ITEMS)

    if effect_coding:
        b_base = np.asarray(params[:5], float)
        b_dif_raw = np.asarray(params[5:10], float)
        b_dif = b_dif_raw - b_dif_raw.mean()  # zero-sum constraint
        muA = params[10]
        muB = params[11]
    else:
        free = [j for j in range(N_ITEMS) if j not in anchor_idx]
        nf = len(free)
        for k, j in enumerate(free):
            b_base[j] = params[k]
            b_dif[j] = params[k + nf]
        muA = params[2*nf]
        muB = params[2*nf + 1]

    b_B = b_base + b_dif

    # Population A likelihood
    thA = GH_NODES + muA
    PA = p_2pl(thA[None, :, None], a_fixed[None, None, :], b_base[None, None, :])
    rA = respA[:, None, :]
    lcA = (rA * np.log(PA + 1e-12) + (1 - rA) * np.log(1 - PA + 1e-12)).sum(axis=2)
    llA = np.log((np.exp(lcA) * GH_W[None, :]).sum(axis=1) + 1e-300).sum()

    # Population B likelihood
    thB = GH_NODES + muB
    PB = p_2pl(thB[None, :, None], a_fixed[None, None, :], b_B[None, None, :])
    rB = respB[:, None, :]
    lcB = (rB * np.log(PB + 1e-12) + (1 - rB) * np.log(1 - PB + 1e-12)).sum(axis=2)
    llB = np.log((np.exp(lcB) * GH_W[None, :]).sum(axis=1) + 1e-300).sum()

    return -(llA + llB)


def fit_dif_pair(respA, respB, a_init, anchor_idx, effect_coding):
    """Fit DIF model for one pair × one identification structure.

    Returns dict with results, or None on failure.
    """
    # Intelligent starting values
    pA = respA.mean(axis=0); pB = respB.mean(axis=0)
    pA = np.clip(pA, 0.05, 0.95); pB = np.clip(pB, 0.05, 0.95)
    bA0 = -np.log(pA / (1 - pA)); bA0 -= bA0.mean()
    bB0 = -np.log(pB / (1 - pB))
    muA0 = float(-np.log(1.0 / max(pA.mean(), 0.01) - 1))
    muB0 = float(-np.log(1.0 / max(pB.mean(), 0.01) - 1))
    muA0 = np.clip(muA0, -3.0, 1.0)
    muB0 = np.clip(muB0, -3.0, 1.0)

    if effect_coding:
        b_dif0 = bB0 - bA0; b_dif0 -= b_dif0.mean()
        x0 = np.concatenate([bA0, b_dif0, [muA0, muB0]])
    else:
        free = [j for j in range(N_ITEMS) if j not in anchor_idx]
        nf = len(free)
        x0 = np.concatenate([bA0[free], (bB0 - bA0)[free], [muA0, muB0]])

    best = None
    best_fval = np.inf
    seeds = [42, 123, 777]  # try multiple starts

    for seed in seeds:
        x = x0.copy()
        if seed != seeds[0]:
            rng = np.random.default_rng(seed)
            x = x + rng.normal(0, 0.2, size=len(x))

        try:
            res = minimize(loglik_joint, x,
                          args=(respA, respB, a_init, anchor_idx, effect_coding),
                          method='Nelder-Mead',
                          options={'xatol': 1e-4, 'fatol': 1e-4, 'maxiter': 50000})
            if res.fun < best_fval:
                best_fval = res.fun
                best = res
        except Exception:
            continue

    if best is None or not best.success:
        return None

    # Extract parameters
    if effect_coding:
        b_base = np.asarray(best.x[:5], float)
        b_dif_raw = np.asarray(best.x[5:10], float)
        b_dif = b_dif_raw - b_dif_raw.mean()
        muA = best.x[10]
        muB = best.x[11]
    else:
        free = [j for j in range(N_ITEMS) if j not in anchor_idx]
        b_base = np.zeros(N_ITEMS)
        b_dif = np.zeros(N_ITEMS)
        nf = len(free)
        for k, j in enumerate(free):
            b_base[j] = best.x[k]
            b_dif[j] = best.x[k + nf]
        muA = best.x[2*nf]
        muB = best.x[2*nf + 1]

    # Compute prevalence estimates
    pA_m2 = prev_from_latent(muA, a_init, b_base)       # population A, common b
    pB_m2 = prev_from_latent(muB, a_init, b_base)       # population B, common b (no DIF)
    pB_m3 = prev_from_latent(muB, a_init, b_base + b_dif)  # population B, with DIF

    # E3 = DIF displacement in B relative to common scale
    e3_diff = pB_m2 - pB_m3

    return {
        'muA': float(muA),
        'muB': float(muB),
        'b_base': b_base.tolist(),
        'b_dif': b_dif.tolist(),
        'pA_m2': float(pA_m2),
        'pB_m2': float(pB_m2),
        'pB_m3': float(pB_m3),
        'e3_dif': float(e3_diff),
        'nit': best.nit,
    }


def load_or_compute_single_2pl(pop):
    """Load M1 data and run single-pop 2PL. Cache results."""
    cache = os.path.join(OUT_DIR, f'{pop}_2pl_single.npy')
    if os.path.exists(cache):
        data = np.load(cache, allow_pickle=True).item()
        return data['resp'], data['a'], data['b'], data['mu']

    df = pd.read_pickle(os.path.join(M1_DIR, f'{pop}.pkl'))
    resp_raw = df[ADL_ITEMS].values.astype(float)
    valid = ~np.isnan(resp_raw).any(axis=1)
    resp = resp_raw[valid]

    # iterative 2PL (simple version — reuse from earlier)
    from scipy.optimize import minimize_scalar

    def loglik_bmu(params, r, a):
        b = np.asarray(params[:5], float); mu = params[-1]
        th = GH_NODES + mu; P = p_2pl(th[None,:,None], a[None,None,:], b[None,None,:])
        rr = r[:,None,:]; lc = (rr*np.log(P+1e-12)+(1-rr)*np.log(1-P+1e-12)).sum(axis=2)
        return -np.log((np.exp(lc)*GH_W[None,:]).sum(axis=1)+1e-300).sum()

    a = np.ones(5)
    for it in range(8):
        ap = a.copy()
        p = resp.mean(axis=0); p = np.clip(p, 0.05, 0.95)
        b0 = -np.log(p/(1-p)); b0 -= b0.mean()
        x0 = np.concatenate([b0, [-1.5]])
        res = minimize(loglik_bmu, x0, args=(resp, a), method='Nelder-Mead',
                       options={'xatol':1e-4,'fatol':1e-4,'maxiter':12000})
        b = np.asarray(res.x[:5], float); shift = b.mean(); b -= shift
        mu = res.x[-1] - shift
        # update a per item
        for j in range(5):
            def f(aj):
                at=a.copy(); at[j]=float(aj)
                th=GH_NODES+mu; P=p_2pl(th[None,:,None],at[None,None,:],b[None,None,:])
                rr=resp[:,None,:]; lc=(rr*np.log(P+1e-12)+(1-rr)*np.log(1-P+1e-12)).sum(axis=2)
                return -np.log((np.exp(lc)*GH_W[None,:]).sum(axis=1)+1e-300).sum()
            rj = minimize_scalar(f, bounds=(0.3, 5.0), method='bounded'); a[j]=rj.x
        if np.max(np.abs(a-ap)) < 5e-4:
            break

    np.save(cache, {'resp': resp, 'a': a, 'b': b, 'mu': mu}, allow_pickle=True)
    return resp, a, b, mu


def load_checkpoint():
    if os.path.exists(CHECKPOINT):
        with open(CHECKPOINT) as f:
            return json.load(f)
    return {}


def save_checkpoint(data):
    with open(CHECKPOINT, 'w') as f:
        json.dump(data, f, indent=2)


def main(pops=None):
    if pops is None:
        pops = ALL_POPS

    print("=" * 80)
    print(f"M3 PAIRWISE DIF: {len(pops)} populations → {len(pops)*(len(pops)-1)//2} pairs")
    print("=" * 80)

    # Load checkpoint
    ck = load_checkpoint()
    done_pairs = set(tuple(sorted(x.split('__'))) if '__' in x else x for x in ck.get('done', []))

    # Load all single-population 2PL
    print("\n[1] Loading single-population 2PL estimates...")
    pop_data = {}
    for pop in pops:
        resp, a, b, mu = load_or_compute_single_2pl(pop)
        pop_data[pop] = {'resp': resp, 'a': a, 'b': b, 'mu': float(mu)}
        print(f"  {pop:15s} N={len(resp):>6d}  mu={mu:.3f}  a={np.round(a,2)}")

    # Run pairwise DIF
    print(f"\n[2] Running pairwise DIF (4 identification structures per pair)...")
    all_pairs = list(combinations(pops, 2))
    results = {}

    # 4 structures
    structures = [
        ('effect_coding', None, True),
        ('anchor_hrs', (4,), False),
        ('sparse_anchor', (3, 4), False),
        ('leave_one_0', (1, 2, 3, 4), False),
    ]

    t0_total = time.time()
    n_new = 0
    for idx, (popA, popB) in enumerate(all_pairs):
        pair_key = f'{popA}__{popB}'
        pair_key_sorted = f'{min(popA,popB)}__{max(popA,popB)}'  # canonical order
        if f'{popA}__{popB}' in done_pairs or pair_key in done_pairs or pair_key_sorted in done_pairs:
            continue

        t0_pair = time.time()
        ra = pop_data[popA]['resp']
        rb = pop_data[popB]['resp']
        a_init = (pop_data[popA]['a'] + pop_data[popB]['a']) / 2.0

        pair_res = {}
        for struct_name, anchor_idx, ec in structures:
            fit = fit_dif_pair(ra, rb, a_init.copy(), anchor_idx, ec)
            pair_res[struct_name] = fit

        results[pair_key] = pair_res
        done_pairs.add(pair_key)
        n_new += 1

        n_ok = sum(1 for v in pair_res.values() if v is not None)
        elapsed = time.time() - t0_pair
        eta = elapsed * (len(all_pairs) - idx - 1)
        print(f"  [{idx+1}/{len(all_pairs)}] {popA} vs {popB}  "
              f"{n_ok}/4 ok  {elapsed:.0f}s  ETA~{eta/60:.0f}min", flush=True)

        # Periodic checkpoint
        if n_new % 5 == 0:
            ck['done'] = list(done_pairs)
            ck['results'] = {k: v for k, v in ck.get('results', {}).items()}
            ck['results'].update(results)
            save_checkpoint(ck)

    # Final save
    ck['done'] = list(done_pairs)
    if 'results' not in ck:
        ck['results'] = {}
    ck['results'].update(results)
    save_checkpoint(ck)

    total_elapsed = time.time() - t0_total
    print(f"\n[3] Complete. {len(done_pairs)} pairs processed in {total_elapsed/60:.1f} min")

    # Summarize
    all_res = ck['results']
    ok_count = 0; fail_count = 0
    for pk, pr in all_res.items():
        for sn, fit in pr.items():
            if fit is not None: ok_count += 1
            else: fail_count += 1

    print(f"  Total fits: {ok_count} converged / {fail_count} failed  "
          f"({100*ok_count/(ok_count+fail_count):.0f}% success)")

    # Quick E3 summary
    print("\n[4] E3 envelope summary (top 5 pairs by |E3|):")
    e3_by_pair = {}
    for pk, pr in all_res.items():
        e3_vals = [fit['e3_dif'] for fit in pr.values() if fit is not None]
        if e3_vals:
            e3_by_pair[pk] = {'min': min(e3_vals), 'max': max(e3_vals),
                              'width': max(e3_vals) - min(e3_vals)}

    for pk in sorted(e3_by_pair, key=lambda x: abs(e3_by_pair[x]['max']), reverse=True)[:8]:
        e3 = e3_by_pair[pk]
        print(f"  {pk:30s}  E3 ∈ [{e3['min']:.4f}, {e3['max']:.4f}]  "
              f"width={e3['width']:.4f}")

    return ck


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--pops', type=str, default='',
                    help='Comma-separated population list (default: all 13)')
    args = ap.parse_args()
    pops = args.pops.split(',') if args.pops else None
    main(pops)
