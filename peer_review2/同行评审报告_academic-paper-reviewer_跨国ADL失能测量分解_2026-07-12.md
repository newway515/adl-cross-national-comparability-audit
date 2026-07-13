# 多视角同行评审报告 —— 跨国 ADL 失能测量分解研究（V8.1）

> 评审引擎：`academic-paper-reviewer` v1.4（full 模式，7 agent / 3 阶段）
> 评审对象：`HANDOFF/` 交接包所述研究协议「How Much of Cross-National ADL Disability Difference Is Measurement Rather Than Health?」（V8.1，pre-results 阶段）
> 评审日期：2026-07-12
> 语言：正文中文，学术术语保留英文
> 证据锚点：`HANDOFF.md`（含 V8.1 定义、M2 卡点、坑清单）、`M1_extract_audit.csv`（真实抽取）、`C8_full_results.csv` / `C8_H2_results.csv`（合成体检）。根目录 `V8.1 决定书`、`C1–C10`、`adl_*.py` 不在评审可见范围，凡涉及其内部按交接自述评估并标注。

---

## 校准声明（重要）

本稿件是**研究协议 / 预结果方案**，主分析结果尚未产出（管线卡在 M2）。因此本评审按「**Registered-Report 式协议评审**」执行：重点评估科学问题的价值、识别策略的合理性、分析方案的严谨性与可行性，而非既成结果。评分标准按目标期刊（Population Health Metrics / BMC Medical Research Methodology / IJE 方法向）的协议/方法学标准校准。

---

# Phase 0 — 领域分析与评审人配置（field_analyst_agent）

## 稿件画像

| 维度 | 判定 |
|---|---|
| 主学科 | 人口健康测量学 / 临床流行病学（population health metrics） |
| 次学科 | 心理测量学（IRT / DIF）、老年流行病学、部分识别方法学 |
| 研究范式 | 观察性二次数据分析 + 测量学分解 + 部分识别 |
| 方法类型 | 6 国际协调老年队列（Gateway Harmonized ADL）、自实现 2PL 边际最大似然、复杂抽样设计方差全传播、四识别结构包络 |
| 目标期刊层级 | 专科方法学期刊（非顶刊综合类） |
| 稿件成熟度 | 协议 / pre-results（M1 完成、C1–C10 审计完成、C8 合成体检通过、M2–M6 未跑） |
| 核心估计量 | E1/E2 点识别；E3/E4 非点识别（包络）；决策量 R_g |

## 五位评审人配置卡

- **EIC（主编）**：*Population Health Metrics* 类期刊主编，全球疾病负担（GBD/WHO）失能测量与跨国可比性专家。关注：期刊契合度、原创性净增量、对负担排名使用者的实际意义、是否过度声称。
- **Reviewer 1（方法学）**：心理测量学 / IRT-DIF 与复杂抽样统计双背景。特别关注：DIF 与 impact 的可识别性、局部独立假设、稀疏二值条目下的参数稳定性、设计方差传播、确证性假设的逻辑自洽。
- **Reviewer 2（领域）**：跨国老年队列失能流行病学专家，熟悉 Gateway to Global Aging harmonized 数据族（HRS/ELSA/SHARE/CHARLS/LASI/MHAS）与既有测量不变性文献。关注：文献占位、Gateway 协调机制、结局定义、领域增量。
- **Reviewer 3（跨学科 / 视角）**：测量不变性与部分识别方法学 + 卫生政策视角（GBD 排名与长期照护规划的下游使用者）。关注：跨学科方法借用（alignment、anchoring vignettes、Manski bounds）、政策可操作性、决策量设计。
- **Devil's Advocate（魔鬼代言人）**：专职挑战核心论点，检测樱桃采摘、确认偏误、逻辑断裂与「So what?」。

> 说明：按技能 checkpoint 规则，Phase 0 配置卡通常交用户确认后再进入 Phase 1。鉴于用户已明确要求一次性产出评审意见文档，此处将配置卡随报告同附，如需调整评审人身份可提出，我据此重跑 Phase 1。

---

# Phase 1 — 五视角独立评审

---

## 评审报告 1／EIC（主编）

### Reviewer Role
EIC — Population Health Metrics 类期刊主编

### Review Focus
期刊契合度、原创性净增量、对负担排名使用者的现实意义、是否存在过度声称。不深入方法细节（留给 R1）。

### Overall Assessment

**Recommendation:** ☒ Major Revision

**Confidence Score:** 4/5

**Summary Assessment：** 本协议要回答一个对 GBD/WHO 失能负担排名有真实政策意义的问题：跨国 ADL 失能率差异中有多少是测量假象而非真实健康。作者用 6 个协调老年队列、13 人群，设计了 M0→M3 四层分解并给出决策量 R_g。作者的科研诚信态度突出——主动做先例碰撞矩阵、拒绝对非识别量声称点估计、拒绝把 E1 浚探成非零。但作为主编我必须指出一个致命张力：作者的核心卖点是「条目/评分/DIF 三层分解」，而其 M1 实证（`M1_extract_audit.csv` + HANDOFF 第 5 节）已表明 E1≈0、E2 仅 CHARLS 非零——**三层里两层在协调数据上几乎空置**，主叙事被自身数据掏空。这不是致命到必须 Reject，因为把「协调已消除 E1/E2」翻转成主发现后，仍有一篇诚实、有用的方法学论文。但当前框架与标题的错配必须先解决。给予 Major Revision。

### Strengths

**S1：政策相关性明确。** 协议第 1 节把研究问题锚定在 GBD/WHO 负担排名的可比性前提上，并给出条件语句（R_g≥1 则该对排序不宜用于负担排名）。对本刊读者群（负担测量、卫生规划）而言，问题本身值得发表。

**S2：科研诚信标准高于同类投稿。** 作者预置 C9 先例碰撞矩阵、点名 Chan 2012 / Luy 2023 / Lee 2018 并**自行禁用「首次/填补空白」表述**（HANDOFF 第 1 节）。这是本刊少见的自律，编辑部乐见。

**S3：可复现性基础扎实。** C10 冻结快照（环境/代码/种子）、C8 前置合成体检、自实现全透明 2PL，均利于复现——契合本刊对 metrics 类论文的透明度要求。

### Weaknesses

**W1：标题与实证图景错配（原创性风险）。**
**问题：** 标题主打「Layered Decomposition of Item-Definition, Scoring, and DIF Effects」，但 HANDOFF 第 5 节已确认 E1≈0、E2 仅 CHARLS。三层实际塌缩为单层 DIF。
**影响：** 审稿人与读者会质疑「为什么三层里两层是空的还要以三层为题」，并怀疑标题夸大。
**建议：** 重构标题与框架，把「协调数据已消除条目/评分层、残余不可比性几乎全为 DIF」作为**主发现**，而非隐藏在方法脚注。
**Severity：** Critical

**W2：决策交付物可能大面积「不可判定」。**
**问题：** R_g 依赖非点识别的 E3/E4，只能以包络报告；叠加 SHARE 小样本，多数国家对的 R_g 区间很可能跨越 1。
**影响：** 论文最终可能对多数国家对给不出明确判定，削弱「决策量」的卖点。
**建议：** 预先声明预期的「不可判定比例」，并改用部分识别区间宽度作为主交付（见 R3）。
**Severity：** Major

**W3：过度声称风险仍存。**
**问题：** 即便自禁「填补空白」，R_g≥1 阈值判定语气仍偏强（暗示可操作的排名弃用规则）。
**建议：** 讨论区显式标注该规则的识别依赖与阈值任意性。
**Severity：** Minor

### Detailed Comments
- **Title & Abstract：** 标题需与 W1 一致重写；建议摘要首句即承认「in harmonized data, item/scoring harmonization is largely complete upstream」。
- **Introduction：** 动机充分，政策锚点清晰；但需把「三层」预期贡献与 M1 实证对齐。
- **Significance：** 对本刊契合度良好，前提是重构后不再夸大三层。

### Questions for Authors
1. 若 E1≈0、E2 仅 CHARLS 是既成事实，论文的**单句核心贡献**改为什么？请用一句话回答。
2. 预计有多大比例的国家对，其 R_g 包络会跨越 1（即不可判定）？这个比例本身是否就是主结果？

### Dimension Scores

| Dimension | Score | Descriptor | Notes |
|---|---|---|---|
| Originality (20%) | 54 | Weak | 三层塌缩为单层 DIF，与既有不变性文献重合上升 |
| Methodological Rigor (25%) | 60 | Adequate | 工程与体检严谨，但核心估计量非识别 |
| Evidence Sufficiency (25%) | 55 | Weak/Adequate | 协议阶段，主结果未出 |
| Argument Coherence (15%) | 58 | Weak | H2 与 M1 结论存在活矛盾 |
| Writing Quality (15%) | 78 | Strong | 交接文档极清晰 |
| **Weighted Average** | **59.5** | **Major Revision** | |

---

## 评审报告 2／Reviewer 1（方法学）

### Reviewer Role
Peer Reviewer 1 — 心理测量学（IRT/DIF）与复杂抽样统计

### Review Focus
DIF 与真实潜特质差异（impact）的可识别性、2PL 局部独立假设、稀疏二值条目下参数稳定性、设计方差传播、确证性假设 H2 的逻辑自洽。

### Overall Assessment

**Recommendation:** ☒ Major Revision（接近 Reject 边界，取决于作者能否接受「核心估计量仅部分可识别」的定位）

**Confidence Score:** 5/5

**Summary Assessment：** 从测量统计角度，本协议在工程与模拟层面做得比多数跨国比较投稿都扎实：自实现 2PL（Gauss-Hermite 21 节点）、C8 用 12 情景×1000 次做前置闸门、设计方差全管线重抽样。作者也正确地承认 E3/E4 非点识别、只报包络。但我必须严肃指出三点会实质威胁「唯一有内容的 DIF 层」可靠性的问题：(1) 5 个二值条目在无假定不变锚的前提下，DIF 与 impact 数学上不可分——这是天花板而非可修的 bug；(2) ADL 条目的局部依赖几乎必然违反 2PL 局部独立，C8 的 S11 已自曝 E4 覆盖率掉到 0.741；(3) 确证假设 H2 在 E1≈0 下必然被自身数据证否，却尚未回收。这些不否定项目价值，但要求把定位从「分解出各层贡献」下调为「在部分识别下界定测量成分的区间」。

### Strengths

**S1：前置模拟闸门是范本级实践。** `C8_full_results.csv` 显示在无失配情景（S1–S9、S12）下模型成功率 100%、CI 覆盖率名义（E3_cover 0.928–0.962）。先在已知真值的合成数据上给算法做体检，再动真实数据，优于绝大多数同类研究。

**S2：对非识别性的诚实处理方向正确。** 承认 E3/E4 只能以四识别结构（effect coding / HRS anchor / sparse robust anchor / leave-one-item）包络报告，并禁用「稳健/点估计」措辞（HANDOFF 第 1 节 V8.1 处理 1），方法学上是对的。

**S3：设计方差全传播。** 采用完整管线重抽样 ≥500 次传播复杂抽样设计方差，而非仅报 iid 标准误，符合 metrics 类论文标准。

### Weaknesses

**W1：DIF/impact 不可识别是硬天花板，需在设计层正视而非留作 limitation。**
**问题：** 5 个二值 ADL 条目，无假定不变锚、无可用 vignettes（C7 已把 SHARE vignettes 判为三重错配降级），则「人群共同阈值平移（DIF）」与「真实潜在健康均值差（impact）」不可同时识别。这正是作者在 V8.1 处理 1 中包络化 E3/E4 的根因。
**影响：** R_g 因此是区间；对多数国家对判定不确定。核心交付物强度受限。
**建议：** 把交付物形式化为 Manski 式部分识别 bounds，报告 identified set 宽度与「是否覆盖 0/覆盖观测差异」，而非阈值判定。
**Severity：** Critical

**W2：局部依赖违反 2PL，直接威胁唯一有内容的 E3 层。**
**问题：** 穿衣/洗澡/如厕在生理上高度成簇，Yen's Q3 残差相关几乎必然显著。C8 的 S11（LD 情景）已显示 E4 覆盖率掉到 0.741、E1 覆盖率 0.761——而这还是在已知真值的合成数据上。
**影响：** 真实数据的局部依赖只会更脏，DIF 估计与包络宽度不可信。
**建议：** 在真实数据上强制做 Q3 / 残差相关诊断；据结果考虑 testlet / bifactor 模型或明确降级报告；把 LD 敏感性纳入主分析而非附录。
**Severity：** Critical

**W3：确证假设 H2 与 M1 结论自相矛盾。**
**问题：** H2「条目定义效应 ≥ DIF 效应」在 E1≈0、E3>0 下必然为假。`C8_H2_results.csv` 显示 H2 在合成数据上恒 holds（因合成数据人为注入大 E1，median_absE1≈3.07–4.72pp），而真实数据图景相反。作者在 HANDOFF 第 5 节承认 E1≈0，却未回收 H2。
**影响：** 一个在数据抽取阶段就已知会被证否的确证假设，会让审稿人立即质疑确证结构的可信度。
**建议：** 撤掉 H2，或重述为「在协调数据上 E1 已被上游消除，故 H2 不再可检验，转为描述性发现」。
**Severity：** Critical

**W4：稀疏条目使题目参数不稳。**
**问题：** `M1_extract_audit.csv` 显示进食困难阳性率极低（CHARLS eata 361/6967≈5.2%、HRS 402/8833≈4.6%）。5 条目中含极稀疏项，2PL 区分度/难度在稀疏格子上估计不稳，DIF 检验功效低。
**建议：** 报告各条目信息函数；对极稀疏条目做敏感性（合并或留一），确认结论不依赖单一稀疏条目。
**Severity：** Major

**W5：SHARE 单国样本对 DIF 检验欠功效。**
**问题：** SHARE 各国 N 1141–2344、阳性数 218–412（`M1_extract_audit.csv`）。5 条目、~300 阳性估单国 DIF，CI 很宽，包络被进一步撑宽。
**建议：** 报告各人群 DIF 估计的功效/CI 宽度；对小样本人群明示判定力受限。
**Severity：** Major

### Detailed Comments
- **Methodology：** 自实现 2PL 需报告收敛诊断、起始值敏感性、Gauss-Hermite 节点数敏感性（21 节点对极端参数是否足够）。
- **Results：** 建议主表并列「点识别层（E1/E2）」与「部分识别层（E3/E4 区间）」，视觉上不要让读者误读包络为点估计。
- **S10 情景警示：** `C8_full_results.csv` 中 S10（2PL 误设）E1_cover=0.0、E4_cover=0.0，说明模型误设下点估计崩溃——真实数据必须报告模型拟合优度与误设敏感性。

### Questions for Authors
1. 在无假定不变锚、无可用 vignettes 的前提下，你如何向审稿人论证 E3 包络的**下界**不是人为收缩得到的？（这正是 H2 抗自证锁想防的，但 H2 已自证否。）
2. 是否已在真实数据上做过 Q3 局部依赖诊断？若 ADL 条目簇显著相关，包络宽度会扩大多少？
3. 自实现 2PL 与任一已发表实现（如 mirt 的等价设定）是否做过交叉验证？

### Dimension Scores

| Dimension | Score | Descriptor | Notes |
|---|---|---|---|
| Originality (20%) | 55 | Weak | DIF 分析与 Chan 2012 高度相邻 |
| Methodological Rigor (25%) | 56 | Weak | 工程严谨但核心非识别 + LD 风险 + H2 矛盾 |
| Evidence Sufficiency (25%) | 54 | Weak | 合成体检充分，真实结果未出 |
| Argument Coherence (15%) | 52 | Weak | H2 自相矛盾未回收 |
| Writing Quality (15%) | 76 | Strong | 方案表述清晰 |
| **Weighted Average** | **57.6** | **Major Revision** | |

---

## 评审报告 3／Reviewer 2（领域）

### Reviewer Role
Peer Reviewer 2 — 跨国老年队列失能流行病学（Gateway harmonized data 族）

### Review Focus
文献占位与增量、Gateway 协调机制、结局与人群定义、领域内既有测量不变性工作的对照。

### Overall Assessment

**Recommendation:** ☒ Major Revision

**Confidence Score:** 4/5

**Summary Assessment：** 作者对协调数据族的领域细节掌握扎实（逐库个人 ID、HRS 宽表拆片、SHARE 内部国码、缺失码均已核验，见 HANDOFF 第 7 节），这在跨国投稿里很难得，M1 抽取结果（患病率 15.1%–28.1%）也合理。但从领域文献看，本协议锁定的 3 篇先例（C9）远不足以覆盖「协调老年队列 + ADL + 测量不变性/DIF」这一相当拥挤的文献带。更根本的问题是**数据基质与研究目标错配**：作者用 Gateway Harmonized ADL（其存在目的就是统一条目与评分）去量化条目定义效应与评分效应——这在设计上自我预置了 E1≈0。领域增量因此比自评「中等」更窄。

### Strengths

**S1：协调数据工程扎实可信。** HANDOFF 第 7 节 9 条坑均已核验落盘（HRS `r14agey_e`、`__c9` 拆片三表连接、SHARE country 内部码 12/28/35/34/29/15/16/25），领域读者可信赖其抽取正确性。

**S2：CHARLS 4 级量表的处理有领域洞见。** 正确识别出 CHARLS 是唯一 4 级量表起点，且在「任何困难=1」二值结局下折叠后患病率与 Harmonized 版相等，故其增量属 E2（评分）而非 E1（条目）——这是精细且正确的领域判断（HANDOFF 第 5 节）。

**S3：人群定义透明。** F2 冻结（目标波次 ≥65、横截面权重 >0、五条目至少一项非缺失）+ F10 硬闸门（未加权 N≥500）清晰可复核。

### Weaknesses

**W1：先例占位审计不足，需扩一圈。**
**问题：** C9 仅锁 Chan 2012 / Luy 2023 / Lee 2018。但 anchoring vignettes / HOPIT（King, Murray, Salomon, Tandon 2004 及其在失能、自评健康的大量应用）、alignment optimization（Asparouhov & Muthén 2014）、协调认知数据 HCAP 的跨国 DIF 方法群，均与本协议方法高度同构。
**影响：** 若不显式对照，审稿人会质疑「为什么不用 alignment」「与 vignettes 文献的净增量是什么」。
**建议：** 补一份二次碰撞矩阵，检索 `harmonized aging cohort + ADL + differential item functioning / measurement invariance` 与 `partial identification + cross-national health comparison`，会议摘要与 supplement 一并查（作者此前有漏检 supplement 先例的教训）。**注：具体文献需经 PubMed/Scholar 复核，本评审不臆造 DOI。**
**Severity：** Major

**W2：数据基质与目标错配。**
**问题：** 用 Harmonized ADL 去量化 E1/E2，等于用「已被上游协调掉条目/评分」的数据研究条目/评分效应，结果 E1≈0 是结构性必然。
**影响：** 三层分解的两层先天空置。
**建议：** 要么坦承并把 E1≈0 作为对 Gateway 协调的审计结论（推荐），要么改用**原始国家问卷**（各国条目数、题干、难度锚不同）作为 E1/E2 的真正基质——后者工作量大但能救活分解。
**Severity：** Critical

**W3：跨波次日历年与期间效应未纳入框架。**
**问题：** 各库取不同波次不同年份（CHARLS W4、HRS W14、SHARE W8=2019-10–2020-03 恰在 COVID 前沿，见 C4）。跨库患病率比较把期间效应/疫情冲击混入「真实健康差异」，而这层混杂不在 item/scoring/DIF 分类内。
**建议：** 在讨论显式承认期间效应；若可能，做同期波次敏感性。
**Severity：** Major

### Detailed Comments
- **Literature：** 需新增测量不变性 / vignettes / alignment 三条文献线，并明确本研究相对每条的净增量边界。
- **Methodology（人群）：** LASI R1 宽格式、无统一 wave（据全局记忆），横截面定位需说明；教育变量 LASI 用非标准 ISCED 边界（C5）需在跨国协变量协调处标注。
- **Discussion：** 建议增加「协调之后仍残留多少不可比」这一对海量协调数据使用者有用的结论。

### Questions for Authors
1. 若目标是量化条目/评分效应，为何不使用原始国家问卷而使用已协调数据？是否评估过用原始表重建 E1/E2 的可行性？
2. 与 alignment optimization 和 anchoring-vignette 文献相比，本研究的**领域净增量**具体是什么？

### Dimension Scores

| Dimension | Score | Descriptor | Notes |
|---|---|---|---|
| Originality (20%) | 52 | Weak | 文献带拥挤，基质错配 |
| Methodological Rigor (25%) | 62 | Adequate | 领域工程扎实 |
| Evidence Sufficiency (25%) | 58 | Weak/Adequate | 抽取可信，主结果未出 |
| Argument Coherence (15%) | 60 | Adequate | 领域论证基本连贯 |
| Writing Quality (15%) | 78 | Strong | 清晰 |
| Literature Integration (opt.) | 55 | Adequate/Weak | 先例覆盖不足 |
| **Weighted Average** | **59.6** | **Major Revision** | |

---

## 评审报告 4／Reviewer 3（跨学科 / 视角）

### Reviewer Role
Peer Reviewer 3 — 测量不变性与部分识别方法学 + 卫生政策视角

### Review Focus
跨学科方法借用（alignment / vignettes / Manski bounds）、决策量 R_g 设计、对下游政策使用者（GBD 排名、长照规划）的实际可操作性、选择性生存等框架外混杂。

### Overall Assessment

**Recommendation:** ☒ Major Revision

**Confidence Score:** 4/5

**Summary Assessment：** 本协议最有前途的部分恰恰是作者没当成主卖点的两点：把 E1≈0 翻转为「协调数据可比性审计」，以及用部分识别区间替代阈值判定。作者已经采纳了「包络报告」这一部分识别的正确直觉，却又退回到 R_g≥1 的点式决策规则——这是方法哲学的内部不一致。从政策使用者角度，一个大面积「不可判定」的 R_g 表用处有限；而一个「各国家对真实差异的识别区间宽度 + 是否覆盖 0」的表反而直接可用。此外，一个重要的临床流行病学混杂——选择性生存——完全未进入框架。

### Strengths

**S1：部分识别直觉正确。** 采纳四识别结构包络，是通向 Manski 式 bounds 的正确一步，跨学科方向对。

**S2：政策问题真实且紧迫。** 把测量可比性与 GBD 负担排名、长照资源规划挂钩，问题对政策使用者有真实意义。

**S3：CHARLS 案例的跨学科示范价值。** 唯一保留严重度梯度的库，可作为「评分规则一旦保留梯度 E2 就非零」的天然对照案例。

### Weaknesses

**W1：R_g 决策量设计偏弱，与自身部分识别哲学冲突。**
**问题：** R_g = 成对测量位移 / 全局中位国家间差异，量纲混合；阈值 1 无标定；在部分识别下退化为区间使阈值判定对多数对失效。
**影响：** 「决策量」名义可操作，实则多数对不可判定。
**建议：** 用 identified set 宽度、区间是否覆盖 0（排序方向是否可识别）、区间是否覆盖全部观测差异（是否可能全为测量）替代 R_g≥1。这更诚实也更可发表（部分识别在失能跨国比较少见）。
**Severity：** Major

**W2：选择性生存混杂未处理。**
**问题：** 高死亡率国家（LASI 印度中位 70）的 65+ 存活者是更强健康筛选的幸存者，其较低失能率部分来自选择偏倚，而非测量、也非当代真实健康。此混杂完全不在 item/scoring/DIF 分类内。
**影响：** 「真实健康差异」被选择偏倚污染，分解无法触及。
**建议：** 讨论显式承认；可结合各国预期寿命/死亡率做定性敏感性。
**Severity：** Major

**W3：跨学科方法借用不足。**
**问题：** 未与 alignment、HOPIT/vignettes、Manski bounds 三条成熟跨学科方法对话。
**建议：** 至少在方法/讨论与三者各对话一段，说明为何选包络路线及其相对优劣。
**Severity：** Minor

### Detailed Comments
- **Estimand：** 建议把主估计量从 R_g 改为 bounds 宽度；R_g 可保留为副图。
- **Policy framing：** 交付语句从「该对不可用于排名」改为「无额外假设下真实差异只能界定在 [a,b]；若区间覆盖 0，排序方向不可识别」。
- **Impact：** 「审计协调数据残余不可比性」对成百上千协调数据使用者有直接价值，应前置为贡献之一。

### Questions for Authors
1. 是否考虑过用 Manski 式 bounds 直接表述估计量，从而与你已采纳的包络哲学完全一致？
2. 选择性生存对跨国患病率差异的可能贡献有多大？是否有队列死亡率数据可做敏感性？

### Dimension Scores

| Dimension | Score | Descriptor | Notes |
|---|---|---|---|
| Originality (20%) | 58 | Weak | 审计 + bounds 重构后可回升 |
| Methodological Rigor (25%) | 60 | Adequate | 包络方向对，R_g 偏弱 |
| Evidence Sufficiency (25%) | 55 | Weak/Adequate | 主结果未出 |
| Argument Coherence (15%) | 60 | Adequate | 哲学内部不一致（包络 vs 阈值） |
| Writing Quality (15%) | 78 | Strong | 清晰 |
| Significance & Impact (opt.) | 66 | Adequate | 政策相关但可能多数不可判定 |
| **Weighted Average** | **59.2** | **Major Revision** | |

---

## 评审报告 5／Devil's Advocate（魔鬼代言人）

### Strongest Counter-Argument（最强反论，约 280 字）

**「这篇研究可能在用一套精致的方法学机器，去量化一个它已经亲手证明不存在的东西。」** 论证如下：项目的整个卖点建立在「跨国失能差异可分解为条目、评分、DIF 三层测量效应」这一前提上。但作者自己的 M1 数据已经给出判决——E1≈0（12/13 人群）、E2 仅 CHARLS。也就是说，在作者选定的数据基质（Gateway Harmonized）上，三层里两层**在数据进入分析之前就注定为零**，因为这套数据的存在目的正是把这两层抹平。于是真正剩下的只有 DIF 一层，而 DIF 层又被作者自己诚实地判定为**非点可识别**、只能报包络。把两件事叠加：论文的确定性结论（E1/E2 点估计）恰好落在「必然≈0」的平凡区域，而论文唯一有实质内容的部分（E3）恰好落在「无法点识别」的不确定区域。这不是执行缺陷，而是**选题—基质—估计量三者的结构性错配**：确定的地方没内容，有内容的地方不确定。若不重构，读者读完最可能的收获是「跨国 ADL 差异里的测量成分，大部分我们无法判定」——一个诚实但接近空手的结论。

### Issue List

| 编号 | 级别 | 维度 | 位置 | 描述 |
|---|---|---|---|---|
| DA-1 | **CRITICAL** | 核心论点 | HANDOFF 第5节 + M1 CSV | 三层分解的两层（E1/E2）被自身数据证明≈0，主叙事空心化 |
| DA-2 | **CRITICAL** | 逻辑一致性 | H2 vs 第5节 | 确证假设 H2 在数据抽取阶段即已知会被证否，仍留在确证结构中 |
| DA-3 | **CRITICAL** | 确定性—内容错配 | V8.1 处理1 | 点识别层（E1/E2）无内容；有内容层（E3）非识别，二者不重叠 |
| DA-4 | MAJOR | 樱桃采摘风险 | 四识别结构包络 | 放弃预注册后，「事后选识别结构」的质疑失去最强防线；H2 抗自证锁又已失效 |
| DA-5 | MAJOR | So-what | R_g 交付 | 若多数国家对 R_g 不可判定，决策量的实际决策价值存疑 |
| DA-6 | MINOR | 确认偏误 | C8 体检 | C8 全为合成数据、H2 恒 holds，易被误读为支持真实世界 H2（作者已警示坑5，但论文正文需防读者误读） |

### Ignored Alternative Explanations / Paths
- **替代解释：** 跨国失能率差异的很大一部分可能既非测量、也非「当代真实健康」，而是**选择性生存 + 期间效应**——两者都在三层框架之外，却可能量级不小。若不排除，"measurement vs health" 的二分本身是伪二分。
- **未走的路：** 用原始国家问卷做 E1/E2、用 Manski bounds 做估计量、把研究重定位为「协调数据可比性审计」——任一条都比现框架更自洽。

### Missing Stakeholder Perspectives
- **协调数据的维护者（Gateway/RAND/SHARE 团队）：** 「E1≈0」其实是对他们协调质量的正面审计证据，这一利益相关方视角完全缺席，却是论文最有价值的受众之一。
- **GBD 负担建模者：** 他们需要的不是「某对不可用」，而是「不确定性区间有多宽」——正指向 bounds 重构。

### Observations（非缺陷）
- 作者的诚信操作（禁用「填补空白」、拒绝浚探 E1、包络化非识别量）值得高度肯定；本 DA 的所有 CRITICAL 均指向**框架定位**而非**学术不端**。项目底子干净，问题在叙事与估计量选择，均可通过重构解决。

---

# Phase 2 — 编辑综合与决定（editorial_synthesizer_agent）

## Editorial Decision

### ☒ Major Revision

（不选 Reject：底层数据/工程/审计资产可信，重构后可成一篇诚实有用的方法学论文。不选 Minor：存在 3 项 DA-CRITICAL，涉及核心叙事、确证逻辑与估计量定位的结构性问题，需实质重构并重新评审。按技能规则，DA 存在 CRITICAL 时决定不得为 Accept。）

## Reviewer Summary

| Reviewer | Role | Recommendation | Confidence |
|---|---|---|---|
| EIC | Population Health Metrics 主编 | Major Revision | 4 |
| Reviewer 1 | 心理测量/IRT-DIF + 复杂抽样 | Major Revision（近 Reject 边界） | 5 |
| Reviewer 2 | 跨国老年失能流行病学 | Major Revision | 4 |
| Reviewer 3 | 部分识别方法学 + 卫生政策 | Major Revision | 4 |
| Devil's Advocate | 核心论点挑战 | 3 项 CRITICAL | — |

**加权均分汇总：** EIC 59.5 / R1 57.6 / R2 59.6 / R3 59.2 → 综合约 **59（50–64 区间 → Major Revision）**，与四位独立评审一致。

## Consensus Analysis

### Points of Agreement

**[CONSENSUS-4]（四位评审 + DA 一致）**
1. **三层分解被自身数据掏空（E1≈0、E2 仅 CHARLS）是当前最致命问题。** EIC-W1、R2-W2、DA-1 均指向此；R1、R3 隐含认同。这是本轮 Major Revision 的首要触发。
2. **核心估计量部分不可识别、R_g 多数对可能不可判定。** R1-W1、R3-W1、EIC-W2、DA-5 一致；建议改用部分识别区间。
3. **科研诚信底子干净，问题在框架定位而非不端。** 五方一致肯定（EIC-S2、R1-S2、R2-S1、R3-S1、DA Observations）。

**[CONSENSUS-3]（3/4 评审一致）**
1. **确证假设 H2 与 M1 结论自相矛盾需回收。** R1-W3、DA-2 明确；EIC 隐含。R2/R3 未直接评述。
2. **先例占位不足，需补 alignment / vignettes / bounds 三条文献线。** R2-W1、R3-W3 明确；EIC 关注原创性时隐含。

### Points of Disagreement

**分歧 1：严重度到 Reject 还是 Major Revision？**
- R1 观点：接近 Reject 边界（核心估计量非识别 + LD + H2 矛盾三重叠加）。
- EIC/R2/R3 观点：Major Revision（重构后可救）。
- 分歧类型：Severity 分歧。
- **编辑裁定：** Major Revision。
- **裁定理由：** R1 的三点均为真，但均可通过「重定位为协调数据可比性审计 + Manski bounds 估计量 + 回收 H2」解决，无一要求推翻已完成的可信资产（M1/C8/设计方差）。按保守可救原则取 Major。

**分歧 2：E1≈0 是缺陷还是卖点？**
- DA/EIC 观点：当前框架下是致命缺陷（主叙事空心）。
- R3/R2 观点：翻转定位后是主发现（对协调数据使用者有价值）。
- 分歧类型：方向分歧。
- **编辑裁定：** 二者统一——「在当前三层叙事下是缺陷，重构为审计叙事后是卖点」。裁定要求作者执行重构（R1'）。

## Decision Rationale（约 260 字）

本决定为 **Major Revision**，四位独立评审加权均分 57.6–59.6 高度一致落于 50–64 区间，且 Devil's Advocate 提出 3 项 CRITICAL，按技能规则排除 Accept。核心依据：作者的数据工程、前置模拟体检、设计方差传播与科研诚信均达到甚至超过目标期刊标准（EIC-S1/S3、R1-S1、R2-S1），但存在一个贯穿全稿的结构性错配——论文点识别的部分（E1/E2）已被自身 M1 数据证明≈0，而唯一有实质内容的部分（E3/DIF）又非点可识别（R1-W1、R3-W1、DA-3）。这使当前「三层分解 + R_g 阈值」的框架既在原创性上与既有不变性文献重合（R2-W1），又在交付物上面临大面积不可判定（EIC-W2）。之所以不 Reject，是因为三条重构路径（审计化定位、bounds 化估计量、回收 H2）均能在不废弃任何已完成可信资产的前提下把项目救成一篇自洽、诚实、对协调数据使用者有用的方法学论文。之所以不 Minor，是因为上述属核心叙事与估计量层面的实质改动，须重新评审。

## Required Revisions（Must Fix）

| # | 修订项 | 来源 | 级别 | 章节 | 预估工作量 |
|---|---|---|---|---|---|
| R1 | 重构核心叙事：把 E1≈0/E2 仅 CHARLS 作为「协调数据可比性审计」主发现，重写标题与框架 | EIC-W1, R2-W2, DA-1/DA-3 | Critical | 标题/引言/框架 | 1–2 周 |
| R2 | 回收确证假设 H2：撤销或重述为描述性发现，消除与 M1 的逻辑矛盾 | R1-W3, DA-2 | Critical | 假设/方法 | 2–3 天 |
| R3 | 估计量重构：以 Manski 式部分识别 bounds（区间宽度 / 是否覆盖 0 / 是否覆盖观测差异）替代或统领 R_g≥1 阈值 | R3-W1, R1-W1, EIC-W2, DA-5 | Critical | 估计量/结果 | 1–2 周 |
| R4 | 局部依赖诊断入主分析：真实数据 Q3/残差相关，据结果处理 testlet 或明确降级；LD 敏感性前置 | R1-W2 | Critical | 方法/结果 | 1 周 |
| R5 | 补二次先例碰撞矩阵：alignment（Asparouhov & Muthén 2014）、HOPIT/vignettes（King et al. 2004）、partial identification 三线，经 PubMed/Scholar 复核 | R2-W1, R3-W3 | Major | 引言/讨论 | 3–5 天 |
| R6 | 稀疏条目与小样本稳健性：条目信息函数 + 极稀疏条目留一/合并敏感性 + SHARE 小样本判定力声明 | R1-W4, R1-W5 | Major | 结果/敏感性 | 3–5 天 |
| R7 | 框架外混杂显式化：期间效应/COVID（SHARE W8）与选择性生存，讨论专段 + 可行的定性敏感性 | R2-W3, R3-W2 | Major | 讨论 | 2–3 天 |
| R8 | 透明性补强：依托 C10 冻结做「解盲前分析冻结」时间戳声明 + 全量代码种子释放，弥补放弃预注册的信誉缺口 | DA-4 | Major | 方法/数据可得性 | 1–2 天 |

### Required Item Details（节选关键项）

**R1：重构为「协调数据可比性审计」**
- 问题：三层叙事的两层被自身数据证≈0（DA-1/DA-3）。
- 要求：新标题与引言以「Gateway 协调已消除条目/评分层，残余跨国不可比性几乎全为 DIF 且部分不可识别」为主线；CHARLS 4 级量表作为 E2 唯一非零的对照案例。
- 验收：审稿人读标题与摘要即能正确预期「E1≈0 是主发现而非隐藏缺陷」。

**R2：回收 H2**
- 问题：H2 在 E1≈0 下必然被证否，却仍在确证结构（DA-2）。
- 要求：撤销 H2 或改述为描述性发现；`C8_H2_results.csv`（合成、恒 holds）在正文明确标注「不构成真实世界 H2 证据」。
- 验收：全稿无「在数据抽取阶段已知会被证否的确证假设」。

**R3：bounds 化估计量**
- 问题：R_g 量纲混合、阈值任意、部分识别下退化为区间（R3-W1）。
- 要求：主交付改为识别区间宽度与覆盖判定；R_g 降为副图并标注识别依赖。
- 验收：政策使用者可从主表读出每对「真实差异被界定在何区间、排序方向是否可识别」。

## Suggested Revisions（Should Fix）

| # | 修订项 | 来源 | 优先级 | 章节 |
|---|---|---|---|---|
| S1 | 自实现 2PL 与已发表实现（如 mirt 等价设定）交叉验证 + 收敛/节点敏感性 | R1 详评 | P2 | 方法/附录 |
| S2 | 主表并列「点识别层」与「部分识别层」，防读者误读包络为点估计 | R1 详评 | P2 | 结果 |
| S3 | 增加「协调后残余不可比性」对协调数据使用者的价值段落（面向 Gateway/GBD 受众） | DA Stakeholder, R3 | P2 | 讨论 |
| S4 | LASI 横截面定位与非标准 ISCED 边界（C5）在协变量协调处标注 | R2 详评 | P3 | 方法 |

## Revision Roadmap

### Priority 1 — 结构性重构（预估 3–4 周）
- [ ] R1：重构核心叙事为协调数据可比性审计（重写标题/引言/框架）
- [ ] R2：回收 H2
- [ ] R3：估计量 bounds 化
- [ ] R4：局部依赖诊断入主分析

### Priority 2 — 内容补充（预估 1.5–2 周）
- [ ] R5：二次先例碰撞矩阵（三条文献线，PubMed/Scholar 复核）
- [ ] R6：稀疏条目 + 小样本稳健性
- [ ] R7：期间效应与选择性生存显式化
- [ ] R8：透明性声明 + 代码种子释放
- [ ] S1–S3：交叉验证、双层并列表、审计受众段落

### Priority 3 — 文字与格式（预估 2–3 天）
- [ ] 统一「非识别量禁用点估计」措辞
- [ ] C8 合成数据在正文全部显式标注「非真实世界估计」
- [ ] S4 协变量协调脚注

### Total Estimated Effort
- **Major Revision 总计：** 约 5–7 周（含真实数据 M2–M6 首轮跑通）

## Revision Deadline
- **建议期限：** 2026-09-06（约 8 周）
- **依据：** Major Revision 常规 6–8 周，本轮含真实数据管线首跑，取上限。
- **延期政策：** 如需延期，请于截止前 1 周告知。

## Closing

我们鼓励作者认真对待各位评审的意见并提交实质性修改稿；修改稿将进入下一轮评审。需要强调：本决定的所有 CRITICAL 均指向**框架定位与估计量选择**，而非数据质量或学术诚信——后两者评审组一致给予高度肯定。作者已完成的 M1 数据层、C1–C10 审计、C8 合成体检、设计方差方案均为可信资产，无需推翻；重构的本质是把这些资产对齐到一个自洽、诚实、且不与自身数据矛盾的目标（协调数据可比性审计 + 部分识别区间）上。若采纳 R1–R3 重构，本项目原创性有望从当前 Weak 回升至 Adequate 偏上，可发表性显著改善。

如目标期刊调整，建议方法学重构后优先投 *BMC Medical Research Methodology*（对部分识别 + 审计型贡献更友好），*Population Health Metrics* 作为并列选项。

---

## 附录：五位评审报告完整版

（见上 Phase 1 各节。DA 报告采用专用格式，含最强反论、CRITICAL/MAJOR/MINOR 分级问题清单、被忽略替代路径、缺失利益相关方视角、非缺陷观察。）

---

### 评审局限声明（诚信要求）

1. 本评审对象为协议/预结果方案，主分析结果未产出；对「结果」维度的评分基于可行性证据而非既成结果。
2. 根目录 `V8.1 决定书`、`C1–C10`、`adl_*.py` 原件不在评审可见范围，四识别结构确切定义、包络计算细节按 HANDOFF 自述评估。
3. R5 所列文献线索（King et al. 2004、Asparouhov & Muthén 2014、Manski partial identification 等）为方向性提示，须经 PubMed/Google Scholar/Semantic Scholar 复核后入稿，本评审未臆造具体 DOI；C9 已锁三篇先例依作者碰撞矩阵采信，未独立复核原文。
4. 本报告由 `academic-paper-reviewer` v1.4 多 agent 框架结构化生成，五视角由单一模型分角色扮演，非五位真人评审；结论供作者自我改进参考，不等同期刊正式评审意见。
