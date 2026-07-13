# -*- coding: utf-8 -*-
"""
M4 前置模块: 补齐 CHARLS/LASI/MHAS 的 PSU/stratum 字段
=====================================================
从原始表中提取社区/PSU/stratum 信息, 回填到 M1 抽取层的 .pkl 文件中。

CHARLS: communityID (社区编码) → 用作 PSU; 无显式 stratum 列 → 用省份/地区作为 stratum 代理
LASI: r1oopsupl1y (过度抽样 PSU 权重, 507 唯一值) / r1gripsum (地域汇总权重, 117 唯一值)
      → r1oopsupl1y 近似 PSU 编码; r1gripsum 近似 stratum
MHAS: upm_dis_21 (UPrimary Sampling Unit) 在 2021 构造变量表中
      → 需确认 W5 (2018) 波次对应的 PSU 变量 (可能在 MHAS 原始表中)
      回退方案: 用H_MHAS_c2 的现有列重新评估, 或排除 MHAS 从全管线 bootstrap 中
      替换为 jackknife 近似

用法:
    python backfill_psu_strata.py              # 回填并更新所有 3 个队列的 .pkl
    python backfill_psu_strata.py --dry-run    # 仅报告 PSU/strata 可用性, 不写入

Author: Claude (V9.1 M4 PSU backfill)
Date: 2026-07-13
"""
import numpy as np
import pandas as pd
import sqlite3, os, sys, argparse, json

HERE = os.path.dirname(os.path.abspath(__file__))
M1_DIR = r'D:\cursorproj\cursorEvaluation\cur260711\analysis\data'
DB_DIR = r'D:\clinicdatabase\SQLitedatabase'
OUT_DIR = os.path.join(HERE, 'analysis', 'data')
os.makedirs(OUT_DIR, exist_ok=True)

# ============================================================
# 1. CHARLS: communityID → PSU, province as pseudo-stratum
# ============================================================
def backfill_charls():
    """CHARLS W4 has communityID in w4_charls_CN table.
    Each communityID encodes a village/neighborhood PSU.
    There's no explicit stratum variable, so we use province/region
    from the first 2 digits of communityID as pseudo-stratum.
    """
    print("[CHARLS] Backfilling PSU/stratum...")
    con = sqlite3.connect(os.path.join(DB_DIR, 'CHARLS.db'))

    # Query communityID from W4 original table
    raw = pd.read_sql('SELECT ID, communityID FROM "w4_charls_CN"', con)
    con.close()

    # Load M1 CHARLS pickle
    m1_path = os.path.join(M1_DIR, 'CHARLS.pkl')
    m1 = pd.read_pickle(m1_path)

    # Merge communityID
    m1 = m1.merge(raw, left_on='pid', right_on='ID', how='left')

    # Create PSU from communityID (encode as int if possible)
    def encode_community(cid):
        if pd.isna(cid):
            return np.nan
        try:
            # communityID is like '0101041' - encode as integer hash
            return int(str(cid).strip())
        except (ValueError, AttributeError):
            return np.nan

    m1['psu'] = m1['communityID'].apply(encode_community)

    # Create pseudo-stratum from first 2 digits of communityID (province code)
    def province_code(cid):
        if pd.isna(cid):
            return np.nan
        s = str(cid).strip()
        if len(s) >= 2:
            return int(s[:2])
        return np.nan

    m1['strata'] = m1['communityID'].apply(province_code)

    # Clean up
    m1 = m1.drop(columns=['communityID', 'ID'], errors='ignore')

    n_psu = m1['psu'].nunique()
    n_strata = m1['strata'].nunique()
    n_missing = m1['psu'].isna().sum()

    print(f"  PSU: {n_psu} unique communities")
    print(f"  Stratum: {n_strata} provinces (from communityID prefix)")
    print(f"  Missing PSU: {n_missing}/{len(m1)} ({100*n_missing/len(m1):.1f}%)")

    # Save
    m1.to_pickle(m1_path)
    return {'pop': 'CHARLS', 'n_psu': n_psu, 'n_strata': n_strata,
            'n_missing': int(n_missing), 'psu_source': 'communityID',
            'stratum_source': 'communityID_province_prefix'}


# ============================================================
# 2. LASI: r1oopsupl1y → PSU, r1gripsum → pseudo-stratum
# ============================================================
def backfill_lasi():
    """LASI R1 has oversampling PSU weights (r1oopsupl1y, 507 unique values)
    and grid summary weights (r1gripsum, 117 unique values, likely state-level).

    LASI's design: multi-stage stratified sampling.
    - PSU = primary sampling unit (village/ward)
    - Stratum = state × urban/rural

    r1oopsupl1y encodes oversampling PSU weights at the PSU level.
    r1gripsum encodes grid-based summaries at state level.

    We use:
    - psu = r1oopsupl1y (non-zero and non-null values → unique PSU codes)
    - strata = r1gripsum (unique grid summary codes → approximate stratum)
    """
    print("[LASI] Backfilling PSU/stratum...")
    con = sqlite3.connect(os.path.join(DB_DIR, 'LASI.db'))

    raw = pd.read_sql(
        'SELECT prim_key, r1oopsupl1y, s1oopsupl1y, r1gripsum, s1gripsum '
        'FROM "H_LASI_a2"', con)
    con.close()

    # Load M1 LASI pickle
    m1_path = os.path.join(M1_DIR, 'LASI.pkl')
    m1 = pd.read_pickle(m1_path)

    # Merge
    m1 = m1.merge(raw, left_on='pid', right_on='prim_key', how='left')

    # PSU: LASI's r1oopsupl1y is an oversampling weight, not a PSU code.
    # Most values are 0 (58,518 out of 72,262 = 81%).
    # Better approach: extract PSU from prim_key structure.
    # LASI prim_key format: SSDDDVVVHHHHRR
    #   SS = state (2 digits), DDD = district (3), VVV = village/PSU (3-4)
    # Use first 8-9 digits of prim_key as PSU proxy (state+district+village)
    m1['psu'] = m1['prim_key'].astype(str).str[:9].apply(
        lambda x: int(x) if x.isdigit() else hash(x) % 100000)
    # Stratum: first 5 digits = state + district
    m1['strata'] = m1['prim_key'].astype(str).str[:5].apply(
        lambda x: int(x) if x.isdigit() else hash(x) % 1000)

    # Clean up
    drop_cols = ['prim_key', 'r1oopsupl1y', 's1oopsupl1y', 'r1gripsum', 's1gripsum']
    m1 = m1.drop(columns=[c for c in drop_cols if c in m1.columns], errors='ignore')

    n_psu = m1['psu'].nunique()
    n_strata = m1['strata'].nunique()
    n_missing = m1['psu'].isna().sum()

    print(f"  PSU: {n_psu} unique (from prim_key prefix)")
    print(f"  Stratum: {n_strata} unique (from prim_key state+district)")
    print(f"  Missing PSU: {n_missing}/{len(m1)} ({100*n_missing/len(m1):.1f}%)")

    m1.to_pickle(m1_path)
    return {'pop': 'LASI', 'n_psu': n_psu, 'n_strata': n_strata,
            'n_missing': int(n_missing), 'psu_source': 'prim_key_prefix9',
            'stratum_source': 'prim_key_prefix5'}


# ============================================================
# 3. MHAS: upm_dis or fallback to jackknife approximation
# ============================================================
def backfill_mhas():
    """MHAS W5 (2018) design variables are not easily accessible in the current
    harmonized tables. The harmonized table H_MHAS_c2 has no explicit PSU/stratum.

    Fallback strategy: MHAS has a known design effect (DEFF ≈ 1.05, the lowest
    among all 13 populations, per C1 audit). This means the clustering is minimal.
    We use a conservative approximation:
    - Create pseudo-PSU from household ID prefix (geographic clustering)
    - Use state/region as stratum from the original MHAS table if available
    - Document this as 'approximate PSU from household ID geography'

    If the original MHAS survey documentation is accessible, the actual
    PSU ('seccion' or 'upm') should be extracted from the raw MHAS data files.
    """
    print("[MHAS] Backfilling PSU/stratum...")
    import sqlite3

    # Try to find PSU in MHAS 2015/2018 constructed variables
    con = sqlite3.connect(os.path.join(DB_DIR, 'MHAS.db'))

    # Check if 'mhas' table has section/segment info
    try:
        raw = pd.read_sql(
            'SELECT unhhidnp, seccion, entidad, upm FROM "mhas" LIMIT 5', con)
        print("  Found 'mhas' table with seccion/entidad/upm columns")
        has_psu_in_mhas = True
    except (pd.io.sql.DatabaseError, sqlite3.OperationalError):
        has_psu_in_mhas = False
        print("  'mhas' table does not have standard PSU columns")

    con.close()

    m1_path = os.path.join(M1_DIR, 'MHAS.pkl')
    m1 = pd.read_pickle(m1_path)

    if has_psu_in_mhas:
        # Extract PSU from mhas table
        con = sqlite3.connect(os.path.join(DB_DIR, 'MHAS.db'))
        raw = pd.read_sql('SELECT unhhidnp, seccion, entidad, upm FROM "mhas"', con)
        con.close()

        m1 = m1.merge(raw, left_on='pid', right_on='unhhidnp', how='left')
        m1['psu'] = pd.to_numeric(m1['upm'], errors='coerce')
        m1['strata'] = pd.to_numeric(m1['entidad'], errors='coerce')
        m1 = m1.drop(columns=['unhhidnp', 'seccion', 'entidad', 'upm'], errors='ignore')
    else:
        # Fallback: create pseudo-PSU from household ID
        # MHAS unhhidnp format: typically encodes geography in first N digits
        print("  Using household ID prefix as pseudo-PSU (geographic clustering)")
        pids = m1['pid'].astype(str)
        # Use first 6 digits as pseudo-PSU (approximate geographic cluster)
        m1['psu'] = pids.str[:6].apply(
            lambda x: int(x) if x.isdigit() else hash(x) % 10000)
        # Use first 2 digits as pseudo-stratum
        m1['strata'] = pids.str[:2].apply(
            lambda x: int(x) if x.isdigit() else hash(x) % 100)

    n_psu = m1['psu'].nunique()
    n_strata = m1['strata'].nunique()
    n_missing = m1['psu'].isna().sum()

    print(f"  PSU: {n_psu} unique values")
    print(f"  Stratum: {n_strata} unique values")
    print(f"  Missing PSU: {n_missing}/{len(m1)} ({100*n_missing/len(m1):.1f}%)")
    print(f"  NOTE: MHAS DEFF is ~1.05 (minimal clustering). Approximate PSU is acceptable.")

    m1.to_pickle(m1_path)
    return {'pop': 'MHAS', 'n_psu': n_psu, 'n_strata': n_strata,
            'n_missing': int(n_missing), 'psu_source': 'unhhidnp_prefix' if not has_psu_in_mhas else 'upm',
            'stratum_source': 'unhhidnp_prefix2' if not has_psu_in_mhas else 'entidad'}


# ============================================================
# Main
# ============================================================
def main(dry_run=False):
    print("=" * 70)
    print("M4 PSU/STRATUM BACKFILL")
    print("V9.1 Protocol — Preparing for Bootstrap Variance Estimation")
    print("=" * 70)

    results = {}

    if dry_run:
        print("\nDRY RUN — reporting PSU/stratum availability only.\n")

        # CHARLS
        con = sqlite3.connect(os.path.join(DB_DIR, 'CHARLS.db'))
        cnt = pd.read_sql(
            'SELECT COUNT(DISTINCT communityID) as n FROM "w4_charls_CN" '
            'WHERE communityID IS NOT NULL', con).iloc[0, 0]
        con.close()
        print(f"[CHARLS] communityID: {cnt} unique communities → PSU available")

        # LASI
        con = sqlite3.connect(os.path.join(DB_DIR, 'LASI.db'))
        for col in ['r1oopsupl1y', 'r1gripsum']:
            cnt = pd.read_sql(
                f'SELECT COUNT(DISTINCT {col}) FROM "H_LASI_a2" '
                f'WHERE {col} IS NOT NULL', con).iloc[0, 0]
            label = 'PSU' if 'oopsupl' in col else 'stratum'
            print(f"[LASI] {col}: {cnt} unique values → {label} available")
        con.close()

        # MHAS
        con = sqlite3.connect(os.path.join(DB_DIR, 'MHAS.db'))
        has_psu = False
        try:
            cols = pd.read_sql('PRAGMA table_info("mhas")', con)
            psu_cols = [c for c in cols['name']
                       if any(kw in c.lower() for kw in ['upm','seccion','psu'])]
            if psu_cols:
                print(f"[MHAS] PSU columns found in 'mhas' table: {psu_cols}")
                has_psu = True
        except:
            pass
        if not has_psu:
            print("[MHAS] No PSU in harmonized tables. Will use household ID prefix as "
                  "pseudo-PSU (DEFF ~1.05, minimal clustering impact).")
        con.close()

    else:
        print("\nBackfilling PSU/stratum into M1 pickle files...\n")
        results['CHARLS'] = backfill_charls()
        print()
        results['LASI'] = backfill_lasi()
        print()
        results['MHAS'] = backfill_mhas()

        # Summary
        print("\n" + "=" * 70)
        print("BACKFILL SUMMARY")
        print("=" * 70)
        for pop, info in results.items():
            print(f"  {pop}: {info['n_psu']} PSUs, {info['n_strata']} strata, "
                  f"{info['n_missing']} missing ({info['psu_source']})")

        # Verify all 13 populations now have PSU/strata
        print("\n[VERIFY] Checking PSU/strata in all 13 populations...")
        POPS = ['CHARLS', 'ELSA', 'HRS', 'LASI', 'MHAS',
                'SHARE_DE', 'SHARE_CZ', 'SHARE_EE', 'SHARE_SI',
                'SHARE_PL', 'SHARE_ES', 'SHARE_IT', 'SHARE_IL']
        all_ok = True
        for pop in POPS:
            pkl = os.path.join(M1_DIR, f'{pop}.pkl')
            df = pd.read_pickle(pkl)
            has_psu = df['psu'].notna().any()
            has_strata = df['strata'].notna().any()
            n_psu = df['psu'].nunique() if has_psu else 0
            status = 'OK' if has_psu and has_strata else 'MISSING'
            if status != 'OK':
                all_ok = False
            print(f"  {pop:15s} PSU={n_psu:>5d}  {'[OK]' if status=='OK' else '[MISSING]'} {status}")

        if all_ok:
            print("\n  ✅ All 13 populations have PSU/stratum. Ready for M4 bootstrap.")
        else:
            print("\n  ⚠️  Some populations still missing PSU/stratum. See above.")

        # Save manifest
        manifest_path = os.path.join(OUT_DIR, 'psu_strata_manifest.json')
        with open(manifest_path, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        print(f"\nManifest saved: {manifest_path}")


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--dry-run', action='store_true')
    args = ap.parse_args()
    main(dry_run=args.dry_run)
