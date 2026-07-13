# 跨国 ADL 失能测量可比性审计 · 最终交接包

> **交接日期：** 2026-07-13
> **方案版本：** V9.1 融合终版
> **课题状态：** M1-M4 全管线完成，论文写作可随时开始
> **GitHub：** https://github.com/newway515/adl-cross-national-comparability-audit
> **原始工作空间：** `D:\cursorproj\cursorEvaluation\cur260711`（代码和数据源）
> **当前工作空间：** `D:\ccproj\ccEvaluation\cccharls07`（方案、评审和分析产出）

---

## 给接手人的 30 秒速览

**这个课题在问什么：** 不同国家报告的老年人 ADL（日常生活活动）失能率差异很大（如中国 28%、斯洛文尼亚 15%）。这些差异中，有多少是老人真实健康的不同，有多少只是因为"问卷条目不同、评分规则不同、文化回答倾向不同"而产生的测量假象？

**我们做了什么：** 用 6 个国际协调老年队列（CHARLS 中国 / ELSA 英格兰 / HRS 美国 / LASI 印度 / MHAS 墨西哥 / SHARE 欧洲 8 国）共 13 个人群、64,644 名 65+ 老年人的数据，系统性地审计了 Gateway to Global Aging Data 协调流程之后跨队 ADL 比较中残留的不可比性。

**核心发现：**
1. Gateway 协调已将条目定义层（E1）统一 —— E1≈0。评分效应（E2）仅 CHARLS 有非零值（16.1pp）。
2. SHARE 内部的 8 个欧洲国家之间 ADL 比较非常可靠（DIF 宽度仅 1.7pp）。但跨队列比较（如 CHARLS vs 欧洲国家）的 DIF 宽度高达 7.1pp。
3. 74% 的国家对排序方向稳定，26% 需要以区间而非点估计解读。
4. 识别不确定性（~3.9pp）远大于抽样不确定性（~0.5pp）——"我们能非常精确地测量 DIF 的点估计，但 DIF 的真值在什么假设下是什么——这个不确定性的量级远远超出了抽样误差。"

**下一步：** 论文写作。所有核心数字已就位，V9.1 方案可直接转写为方法节。

---

## 文件夹结构

```
HANDOFF_FINAL/
├── README.md                           ← 本文件：交接包索引入口
├── 下一步做什么.md                      ← 简洁的行动清单
│
├── 01_课题方案/                         ← 方案层的完整演进史
│   ├── 研究方案_V9.1_融合终版.md         ← 【主文件】当前唯一定稿方案（约 60KB）
│   ├── 研究方案_V9.0_协调之后的残差.md    ← 我方初版 V9.0
│   └── 课题研究报告_V9.0-B_同行初版.md    ← 同行独立完成的 V9.0 对比版
│
├── 02_评审与演变/                       ← 评审驱动的课题演进全链条
│   ├── HANDOFF_V8.1_原始交接.md         ← 原始 V8.1 三层分解框架
│   ├── 评审1_叙述版评审.md              ← 第一份独立评审
│   ├── 评审2_ScholarEval评审.md         ← 第二份八维度评分
│   ├── 评审3_临床流病专家评审.md        ← 第三份深度方法学评审
│   ├── 评审4_多视角同行评审.md          ← 第四份五视角多 agent 评审
│   ├── 两份V9.0方案对比分析.md          ← 两份独立 V9.0 的异同分析
│   ├── 附录_评审共识与重构决策对照表.md  ← 四份评审 → V9.1 的逐条映射
│   └── Chan2012_全文逐行对照报告.md     ← 与直接先例的七维度对比
│
├── 03_分析代码/                         ← 所有分析的完整代码
│   ├── smoke_charls_2pl.py              ← CHARLS 冒烟测试（2PL 前置验证）
│   ├── adl_operationalizations.py       ← M2 操作化（CHARLS 4级提取 + 全人群 2PL）
│   ├── run_m3_pairwise.py               ← M3 成对 DIF 分解（78对 × 4结构）
│   ├── finalize_V9_deliverables.py      ← 识别边界图 + 可比性标记 + 论文数据
│   ├── adl_variance_delta.py            ← M4 Delta 方法方差传播
│   ├── backfill_psu_strata.py           ← PSU/stratum 回填到 M1 数据
│   └── adl_variance.py                  ← M4 完整 bootstrap（被 Delta 替代）
│
├── 04_分析结果/                         ← 所有实证产出的最终数据
│   ├── 分析结果总览.md                   ← 【主文件】所有核心数字的汇总表
│   ├── V9_paper_frontier_data.json      ← 论文就绪：78对识别边界 + 全景统计
│   ├── m3_pairwise_checkpoint.json      ← M3 全部 78对 × 4结构 DIF 结果
│   ├── m4_bootstrap_ci_summary.json     ← M4 真 bootstrap CI 总结
│   ├── m2_operationalization_results.json ← M2 操作化结果（含 CHARLS E2）
│   └── m4_bootstrap_2pl_checkpoint.pkl  ← M4 12人群 × 50 reps 全量 bootstrap
│
└── 05_论文材料/                         ← 论文写作的直接输入材料
    ├── 论文写作清单.md                   ← 已完成 / 待完成的对照清单
    ├── 论文核心数字卡片.md               ← 结果节、讨论节所需的所有精确数字
    └── 论文建议结构.md                   ← 推荐的目标期刊、标题、摘要和章节大纲
```

---

## 如何从零了解这个课题（阅读顺序）

| # | 文件 | 说明 |
|---|---|---|
| 1 | `01_课题方案/研究方案_V9.1_融合终版.md` | 完整的研究协议——背景、方法、数据、预期结果、局限，一切都在这里 |
| 2 | `04_分析结果/分析结果总览.md` | 所有核心数字的汇总表——看这一份就知道我们发现了什么 |
| 3 | `下一步做什么.md` | 如果只想继续推进课题，从这里开始 |
| 4 | `02_评审与演变/附录_评审共识与重构决策对照表.md` | 理解为什么方案是现在这个样子 |
| 5 | `05_论文材料/论文核心数字卡片.md` | 写论文时需要的所有精确数字 |

---

## 关键事实速查

| 事实 | 数值 |
|---|---|
| 原始工作空间（代码 + 数据） | `D:\cursorproj\cursorEvaluation\cur260711` |
| M1 抽取数据（13 人群 .pkl） | `D:\cursorproj\cursorEvaluation\cur260711\analysis\data\` |
| SQLite 数据库 | `D:\clinicdatabase\SQLitedatabase\{队列名}.db` |
| Python 环境 | `D:\miniconda\envs\myenv\python.exe` (Python 3.13.11, 无第三方 IRT 包) |
| GitHub 仓库 | https://github.com/newway515/adl-cross-national-comparability-audit |
| 分析产出目录 | `D:\ccproj\ccEvaluation\cccharls07\analysis\data\` |

---

## 课题一句话

**世界协调老年队列中，跨国 ADL 失能差异的可靠性因比较类型而异——SHARE 内部高度可靠，跨队列需以识别区间解读——我们提供了一张识别边界图，告诉使用者具体哪些比较可信。**

---

*交接包完整。下次打开新工作空间时，从本文件夹的 `README.md` 和 `下一步做什么.md` 开始。*
