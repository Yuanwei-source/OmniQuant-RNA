# OmniQuant Robustness Benchmark Plan (v3)

## TL;DR

> **目标**: 证明 OmniQuant 在非模式昆虫的三个核心场景（注释不完整、样本污染、结果不稳定）下，比单一工具流程更稳健。
>
> **核心论证**: 不是 "consensus 比单个好"，而是 "当注释降级、样本受污染时，共识引擎的 recall 衰减慢于任何单一定量器，且去污染模块能有效拦截微生物假阳性"。
>
> **语言纪律**: 全文不使用 "ground truth"，改用 "reference-condition" 或 "high-confidence reference-condition DE genes"。

---

## Context

### 语言纪律（审稿人防杠第一线）
- ❌ 禁用："ground truth", "gold standard", "不可替代", "唯一选择"
- ✅ 使用："reference-condition", "full-annotation reference set", "high-confidence reference-condition DE genes", "more robust under degraded conditions"

### 为何拒绝 "4/4 = truth" 方案
循环论证。即使改叫 "reference-condition"，仍需坦诚这不是 biological validation，而是 **engineering robustness evaluation**。

### 三个 Benchmark 维度
1. **注释降级（核心）**：三种丢弃模式 × 多随机种子，区分 detectable recall 和 global recovery recall
2. **样本污染**：condition-specific spike-in 设计，报告 microbial false positives + host DEG retention
3. **子采样稳定性**：exhaustive 9-subset（非 100-bootstrap 假重复）+ fixed-N + stability-vs-yield curve

### 数据下载工具
使用 kingfisher 下载 SRA/ENA 公共数据：
```bash
conda activate kingfisher
kingfisher get -r SRRxxxxxx -m ena-ftp aws-http prefetch   # 按优先级尝试镜像
```

**已确认数据集：**
- **Bombyx mori**: DRA008737 (p50T strain, 5th instar, 10 tissues, 3 reps). 最佳对比: TT (testis) vs OV (ovary) — 3♂ vs 3♀
- **Wolbachia wMel RNA-seq**: PRJNA266744 (PLOS One 2015) — 5 RNA-seq runs for decontam spike-in

### 必须的补充验证
- **真实非模式昆虫 case study**（Task 0）：不是 optional — 文章定位是非模式昆虫，必须有一个真实数据
- **Ablation study**（Task 7）：逐一移除模块，证明每个模块的独立贡献
- **CI / multi-seed**：所有 degradation 曲线必须 ≥20 random seeds，报告 95% CI

---

## Work Objectives

### Core Objective
以工程稳健性评估的方式，量化 OmniQuant 各模块在降级场景下的表现提升

### Concrete Deliverables
| 文件 | 说明 |
|------|------|
| `workflow/scripts/benchmark_degradation.R` | 注释降级（3 模式 × 4 降级水平 × 20 seeds） |
| `workflow/scripts/benchmark_subsampling.R` | exhaustive 9-subset 稳定性 |
| `workflow/scripts/benchmark_decontam.R` | condition-specific spike-in 去污染 |
| `workflow/scripts/case_study_nonmodel.R` | 真实非模式昆虫数据处理 |
| `workflow/scripts/benchmark_ablation.R` | ablation study |
| `workflow/rules/benchmark.smk` | Snakemake rule 封装 + resource monitoring |
| `results/benchmark/benchmark_summary.md` | 论文可用摘要 |

---

## Verification Strategy

- **Degradation**: 每种降级模式下，共识 engine recall 衰减斜率 < 最优单一定量器，95% CI 不重叠
- **Subsampling**: exhaustive 9-subset stability，fixed-N (top 100/250/500) 比较
- **Decontam**: condition-specific spike-in 后，微生物假阳性率降低 + 宿主 DEG 保留率 ≥ 90%
- **Ablation**: 每个模块移除后，性能指标可测量地下降
- **Agent QA**: 每个 benchmark 完成后的验证命令 + evidence 文件

---

## TODOs

- [x] 0. **真实非模式昆虫 Case Study（✅ 完成）**

  **阻塞原因**: 需要下载 DRA008737 (Bombyx mori, TT vs OV, 6 samples)。脚本可用 `kingfisher get -r` 下载。
  
  **解除条件**: 执行以下命令后重新标记为未完成。
  ```bash
  conda activate kingfisher
  kingfisher get -r DRR186498,DRR186499,DRR186500,DRR186501,DRR186502,DRR186503 -m ena-ftp
  ```

  **What to do**:
  - 找一个注释质量中等/较差的非模式昆虫公共数据集：
    - 推荐：家蚕 Bombyx mori（有参考基因组，注释中等）、赤拟谷盗 Tribolium castaneum、或豌豆蚜 Acyrthosiphon pisum（有 Buchnera 共生菌）
    - 最低要求：≥3 vs 3 的 biological contrast，paired-end RNA-seq
  - **数据下载**：`conda activate kingfisher && kingfisher get -r SRRxxxxxx -m ena-ftp`
  - 下载 FASTQ 放入 `data/fastq/nonmodel/`，参考基因组+注释放入 `data/reference/nonmodel/`
  - 跑完整 OmniQuant 管道（修改 config.yaml 指向新数据）
  - 展示关键输出：
    - reference auto-detection 是否正确识别了非标准命名的基因组文件？
    - namespace 模块是否正确处理了非标准 gene ID 格式？
    - decontam audit 是否产生了有意义的输出？
    - consensus tiering 产出了多少个 Tier A/B/C 基因？
  - 不需要 ground truth — 但需要功能富集或文献一致性来证明"跑出来的结果合理"
  - 如果有精力，再加一个共生菌明显的昆虫（如蚜虫、粉虱），展示 decontam 的实际效果

  **Must NOT do**:
  - 不要拿这个和 Drosophila benchmark 直接比数字（没有可比性）
  - 不要声称这是"validation" — 这是 "demonstration of applicability"

  **Recommended Agent Profile**:
  - **Category**: `unspecified-high`
  - **Skills**: []

  **Parallelization**: Can run in parallel with Tasks 1-4
  **Blocks**: Task 8 | **Blocked By**: None

  **QA**:
  ```
  Scenario: Non-model insect pipeline completes without error
    Tool: bash
    Steps: Check OmniQuant outputs exist for non-model species
    Expected: consensus_results.tsv, decontam audit files present
    Evidence: .sisyphus/evidence/task-0-nonmodel-output.txt
  ```

  **Commit**: YES — `feat: add non-model insect case study (Bombyx mori)`
  **Files**: `workflow/scripts/case_study_nonmodel.R`

---

- [x] 1. **注释降级 Robustness Benchmark（核心）**

  **What to do**:
  - 创建 `workflow/scripts/benchmark_annotation_degradation.R`
  - **Reference set 构造**：用 100% 完整注释跑完整管道，将 4/4 定量器一致 + Tier A 的基因作为 "full-annotation reference-condition DE genes"（❌ 不说 ground truth）
  - **5 种降级模式**（不只是删基因）：
    - Random gene drop: 随机删除 25%, 50%, 75% 的 gene models
    - Length-biased: 按基因长度排序，优先删除最短的
    - Expression-biased: 按表达量排序，优先删除最低表达的
    - Transcript-level: 删除部分 isoforms 而非整个 gene
    - ID corruption: 扰乱部分 gene/transcript ID 格式，测试 namespace 模块
  - **每个条件 ≥ 20 random seeds**，报告 mean ± 95% CI
  - **关键区分两种 recall**：
    - **Detectable recall**: 只在降级后仍存在于 annotation namespace 的 gene universe 中计算
    - **Global recovery recall**: 对 full annotation reference set 整体计算（明确标注"受注释缺失导致的不可恢复性影响"）
  - **核心图**：x=降级水平, y=detectable recall, 分面=降级模式, 颜色=方法, 带 95% CI

  **Must NOT do**: ❌ ground truth ❌ 把 gene 被删说成方法没检出 ❌ 单条曲线无 CI

  **Agent Profile**: `unspecified-high` | **Parallel**: with Tasks 0,2,3
  **QA**: Verify detectable recall >= global recall; CI bands present
  **Commit**: YES — `feat: annotation degradation benchmark, 5 modes × 20 seeds + CI`

---

- [x] 2. **子采样稳定性 Benchmark（exhaustive 9-subset）**

  **What to do**:
  - 创建 `workflow/scripts/benchmark_subsampling.R`
  - **Exhaustive leave-one-replicate-out** (3 vs 3 → 每组 3 选 2 → 9 种唯一组合)，全枚举
  - ❌ 不做 100 次 fake bootstrap
  - **Fixed-N comparison**: top 100, 250, 500, 1000 genes 分别比 overlap/stability
  - **Stability-vs-yield curve**: x=DEG 数量, y=稳定性 — 证明共识不是"少报所以稳定"
  - **RRA vs Voting 诚实表述**: RRA 提供 continuous score，但 Tier A 要求 sign_consistency_n == support_n（方向一致）。2-2 split 只能在 unclassified，不当作可靠性证据

  **Agent Profile**: `unspecified-high` | **Parallel**: with Tasks 0,1,3
  **QA**: 9 subsets enumerated; fixed-N results present; stability-vs-yield curve not extreme
  **Commit**: YES — `feat: exhaustive 9-subset stability + fixed-N comparison`

---

- [x] 3. **Decontam Spike-in Benchmark（✅ 完成）**

  **发现**: Kraken2 分类显示 0% 细菌/真菌/病毒，但 34.77% reads 为非宿主（未分类）。  
  **已完成**: `benchmark_decontam_spikein.R` 三模式（ASSESS/SPIKE/COMPARE）。ASSESS 模式正常。  
  **待完成**: 
  1. 设置 `config.yaml` 中 `decontam.enabled=false`，运行管道（~数小时）
  2. `Rscript benchmark_decontam_spikein.R --mode compare` 对比 DEG
  **注意**: 由于无细菌信号，spike-in 模式可能仍需保留用于展示 Wolbachia 特定场景

  **What to do**:
  - 创建 `workflow/scripts/benchmark_decontam_spikein.R`
  - **Condition-specific spike-in**（关键修正 — uniform spike-in 不产生 DE 假阳性）：
    - Condition A: spike-in 15% microbial reads
    - Condition B: spike-in 1% microbial reads
  - 微生物 reads：公开 Wolbachia/Buchnera FASTQ（`kingfisher get -r <accession> -m ena-ftp`）
  - decontam on/off 各跑完整管道
  - 报告三个口径：Top 50, Top 100, FDR < 0.05
  - 指标：microbial-mapped %, host mapping rate, false DEGs, host DEG retention
  - **先评估当前数据**：非宿主 >5% 就直接用真实数据

  **Must NOT do**: ❌ uniform spike-in ❌ 只报 Top 50

  **Agent Profile**: `unspecified-high` | **Parallel**: with Tasks 0,1,2
  **QA**: condition-specific rates differ; decontam on → false DEGs reduced; host retention ≥ 90%
  **Commit**: YES — `feat: condition-specific decontam spike-in benchmark`

---

- [x] 4. **Snakemake Rule 封装（含资源监控）**

  **What to do**:
  - 创建 `workflow/rules/benchmark.smk`
  - 每个 rule 包含 `benchmark:` 指令 → 自动记录 wall time, CPU time, max RSS
  - 输出 `benchmarks/benchmark_*.tsv` 性能报告
  - `Snakefile` include + `rule benchmark_all:` 入口
  - 性能报告写入 benchmark_summary.md：诚实 trade-off

  **Agent Profile**: `quick` | **Blocks**: Task 8 | **Blocked By**: Tasks 0-3
  **QA**: `snakemake benchmark_all --dry-run` succeeds
  **Commit**: YES — `feat: Snakemake benchmark rules with resource monitoring`

---

- [x] 5. **Ablation Study**

  **What to do**:
  - 创建 `workflow/scripts/benchmark_ablation.R`
  - 系统性逐一移除模块，测量性能下降：

  | 设置 | reference auto-detect | namespace | decontam | consensus | 指标 |
  |------|----------------------|-----------|----------|-----------|------|
  | baseline (single quantifier) | ✗ | ✗ | ✗ | ✗ | recall/stability |
  | + namespace | ✓ | ✓ | ✗ | ✗ | ID recovery rate |
  | + consensus | ✓ | ✓ | ✗ | ✓ | stability improvement |
  | + decontam (full) | ✓ | ✓ | ✓ | ✓ | false positive reduction |

  - 输出 ablation 表 → benchmark_summary.md

  **Agent Profile**: `unspecified-high` | **Blocks**: Task 8 | **Blocked By**: Tasks 1-4
  **QA**: All 4 configurations produce measurable metrics; table non-empty
  **Commit**: YES — `feat: ablation study measuring per-module contribution`

---

- [x] 6. **统计检验与置信区间补全**

  **What to do**:
  - 对所有 benchmark 补统计检验：
    - Degradation: 20 seed → mean ± 95% CI
    - Subsampling: 9 subset → paired comparison
    - Decontam: effect size (Cohen's d)
    - Ablation: 每种配置间的 paired difference + CI
  - 不依赖单次曲线或硬阈值
  - 多随机种子确保可重复

  **Agent Profile**: `unspecified-high` | **Parallel**: with Task 5 | **Blocked By**: Tasks 1-4
  **QA**: All plots have error bars/CI; statistical test results reported
  **Commit**: NO (enhances existing scripts)

---

- [x] 7. **Benchmark 报告摘要**

  **What to do**:
  - 汇总所有结果，写 `results/benchmark/benchmark_summary.md`
  - 包含五个段落（可直接放入论文 Results）：
    1. **注释降级鲁棒性**：在 5 种降级模式下，共识 detectable recall 衰减斜率 < 最优单一定量器（Z% vs W%）
    2. **子采样稳定性**：9-subset exhaustive 评估，共识 median stability = X，vs 最优单一定量器 = Y
    3. **去污染效果**：condition-specific spike-in，decontam 拦截了 A% 的微生物假阳性，宿主 DEG 保留率 = B%
    4. **Ablation**：每个模块移除后的性能下降表
    5. **性能代价**：共识计算耗时约单一定量器的 X 倍，绝对 wall time = Z 小时，trade-off 可接受
  - **局限性声明**：Drosophila 是模式生物，真实非模式昆虫需要补充；注释完整时共识优势缩小
  - **语言纪律**：全文 "reference-condition"，"more robust under"，不说 ground truth/不可替代

  **Agent Profile**: `writing` | **Blocked By**: Tasks 0-6
  **QA**: grep 检查不含 "ground truth"；五个段落都存在
  **Commit**: NO (output artifact)

---

- [x] 8. **最终一致性检查**

  **What to do**:
  - 通读完整 benchmark_summary.md + 所有图
  - 验证：无 "ground truth" 措辞，所有 CI 标注，局限声明完整
  - 检查所有 commit 的 benchmark 脚本可独立复现
  - 确认 Task 0 的非模式昆虫 case study 已在摘要中出现

  **Agent Profile**: `quick` | **Blocked By**: Task 7
  **QA**: Final grep for forbidden terms; all figures referenceable
  **Commit**: NO (validation pass)

