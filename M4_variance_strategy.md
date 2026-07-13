# M4 Bootstrap 方差估计 · 策略调整与技术报告

> **日期：** 2026-07-13
> **冒烟结果：** 3 人群 × 3 对 × 1 次复制 ≈ 344 秒（单核，无并行）
> **外推：** 13 人群 × 78 对 × 1000 次复制 ≈ 344/3 × 78 × 1000 ≈ 8,900,000 秒 ≈ **约 100 天**（单核串行不可行）

---

## 为什么单核 Bootstrap 不可行

每 bootstrap 复制内要重做：
1. 13 个单人群 2PL 拟合（每个 10 次迭代 × 5 条目 profile likelihood）
2. 78 对 pairwise DIF（每对 Nelder-Mead 50,000 次迭代上限）

即使算上 12 核并行，100 天/12 ≈ 8 天仍然太慢。

## 实际可行的策略

**核心思路：不降低方法学严谨性，但缩小 bootstrap 的覆盖范围。**

### 方案 A（推荐）：分层 Bootstrap — 仅对已收敛的 M3 点估计传播方差

不对**全管线**做 bootstrap（M3 pairwise DIF 单次就 7 分钟），而是：

1. **Step 1 — Bootstrap 仅传播单人群 2PL 参数的不确定性**（每人群 ~30 秒/rep）
   - 对 13 人群各做 1000 次 bootstrap 2PL（每 rep ~30 秒，总计 ~390 分钟/12 核 ≈ 33 分钟）
   - 产出：a, b, mu 的 bootstrap 分布 + 患病率估计的 bootstrap 分布

2. **Step 2 — Delta 方法传播到 E3**
   - 使用 M3 pairwise checkpoint 中已收敛的点估计作为锚
   - 将 bootstrap 产生的 2PL 参数不确定性通过 delta 方法（数值雅可比）传播到 DIF 包络
   - **不重跑 pairwise DIF** —— 仅在已收敛的点估计周围做局部线性近似

3. **Step 3 — 报告**
   - 对每个国家对，E3 的 bootstrap CI = 95% 区间（跨 bootstrap rep 的分位数）
   - 识别区间宽度的 bootstrap CI
   - "排序方向不可识别"标签的 bootstrap 稳定性

### 方案 B（备选）：仅 bootstrap 关键对

仅对 policy-relevant pairs 做全管线 bootstrap（如 CHARLS-HRS, CHARLS-ELSA, ELSA-HRS, LASI-HRS 等 ~10 对），其余对用 delta 方法近似。可大幅缩短计算时间，但无法给出全 78 对的 CI。

---

## 对论文的影响

### 如果选择方案 A（推荐）

**方法节写为：** "≥1000 次 bootstrap 重抽样传播复杂调查设计方差至单人群 2PL 参数估计；DIF 包络的方差通过 delta 方法（数值雅可比）从 2PL 参数协方差矩阵传播至患病率量纲的 E3/E4 位移。"

**方法学严谨性：** Delta 方法在已收敛的 MLE 附近是渐近有效的。Bootstrap 传播了参数不确定性（包括复杂调查设计导致的方差膨胀），delta 方法将其转化为 E3 的 CI——比纯 iid SE 更严谨，比全管线全 bootstrap 更可行。

### 如果选择方案 B

可行但论文结果表会有"完整 CI"和"近似 CI"两类——需要在脚注中区分。

---

## 推荐行动

**立即采用方案 A（Delta 方法），当天完成，提交 GitHub。** 之后 M4 CI 就可以写进论文。

具体步骤：
1. 对全 13 人群 × 1000 reps 做 2PL bootstrap（单人群，快速，~30 秒/rep × 13 pop = 6.5 小时后并行）
2. 对每 bootstrap rep 计算 E1/E2/E3/E4（使用已收敛的 M3 参数 + delta 传播）
3. 输出每个国家对的 E3 CI 表
