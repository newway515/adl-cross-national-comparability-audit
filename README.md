# ADL Cross-National Comparability Audit · Frozen Study Protocol

> **Repository created:** 2026-07-12
> **Protocol version:** V9.1 (frozen, pre-results)
> **Status:** All analysis decisions frozen before contact with real-data layer M2+

---

## What this repository contains

This is the **immutable, timestamped, pre-results** study protocol for:

**"Auditing the Comparability of Cross-National ADL Disability Estimates: A Partial-Identification Decomposition of Item, Scoring, and Differential-Item-Functioning Effects, with an Identifiability-Frontier Map"**

The protocol freezes all analysis decisions (F1–F17, documented in Protocol Appendix C) before the main analysis pipeline (M2–M6) is executed on real data.

## Repository structure

```
.
├── README.md                                     ← This file
├── protocol/
│   └── 研究方案_V9.1_融合终版.md                   ← Complete frozen protocol (Chinese)
├── handoff/
│   ├── HANDOFF.md                                ← Original V8.1 handoff document
│   ├── M1_extract_audit.csv                      ← Real-data extraction audit (13 populations, completed)
│   ├── C8_full_results.csv                       ← Synthetic-data simulation results (algorithm validation)
│   └── C8_H2_results.csv                         ← H2 sensitivity simulation results (synthetic)
├── reviews/
│   ├── 课题评审报告_跨国ADL失能测量分解.md          ← Narrative review #1
│   ├── 课题ScholarEval结构化评审报告.md             ← ScholarEval 8-dimension review #2
│   ├── 课题评价_跨国ADL失能测量分解_2026-07-12.md    ← Clinical epidemiology review #3
│   └── 同行评审报告_academic-paper-reviewer.md       ← Five-perspective structured review #4
└── supplementary/
    ├── 两份V9.0方案对比分析.md                      ← Comparison of two independent V9.0 drafts
    ├── 附录_评审共识与重构决策对照表.md              ← Review-to-revision mapping table
    └── README_V9.md                               ← Original V9.0 package index
```

## Key facts about this protocol

1. **Data:** 6 harmonized aging cohorts (CHARLS/ELSA/HRS/LASI/MHAS/SHARE), 13 populations, N=64,644
2. **Core empirical constraint:** In Harmonized data, E1 (item-definition effect) ≈ 0 globally, E2 (scoring effect) ≠ 0 only in CHARLS — the substantive measurement non-comparability residual is entirely E3 (DIF)
3. **E3 is NOT point-identifiable** with 5 binary items and no anchoring vignettes — reported as partial-identification intervals (Manski bounds) across 4 identification structures
4. **All analysis decisions (F1–F17) frozen before M2 execution** — see Protocol Appendix C

## Transparency statement

All analysis decisions documented in this repository were frozen before contact with the real-data analysis layer (M2 onwards). The commit history of this repository, combined with the Zenodo DOI timestamp, provides third-party verifiable evidence that these decisions were not made post-hoc.

The C1–C10 pre-analysis audit reports (design effects, item mapping, education comparability, vignette feasibility, freeze snapshot, etc.) are available as supplementary material upon request — they were completed before this repository was created but require desensitization before public release.

## Citation

If you use or reference this protocol, please cite the Zenodo DOI associated with this repository's v1.0.0 release.

## License

This protocol is made available for transparency and reproducibility purposes. Analysis code will be released under MIT License with the final publication. Raw data are not redistributable — access them via the respective cohort data portals (Gateway to Global Aging Data, CHARLS, LASI, MHAS).

---

*Repository prepared for Zenodo archival. Contact the study team for access to the C1–C10 pre-analysis audit reports and the C10 freeze snapshot (environment/code/seeds).*
