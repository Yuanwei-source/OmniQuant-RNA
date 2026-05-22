# Epicauta impressicornis 滞育转录组 — 部署与执行计划 v2

## TL;DR

> **目标**: 在扁角豆芫菁（Epicauta impressicornis）滞育转录组数据（non-diapause vs diapause, 3 vs 3）上运行 OmniQuant-RNA 两遍（decontam ON + OFF），产出共识 DEA 结果和去污染对比数据，填入论文 §3.6。
>
> **数据源**: NAS `/home/dell/cxslab_nas/yuanwei/diapause/` — 8 个样本已下载（16 个 FASTQ），取其中 6 个核心样本（ND+ D），跳过 DT 组。
>
> **参考基因组**: `data/reference/epicauta/genome.fasta` + `annotation.gff3` ✅ 已就位
>
> **预计耗时**: 符号链接 ~1min + 管道运行 ~8h × 2 = ~16h

---

## Context

### 数据状态

| 项目 | 状态 |
|------|------|
| FASTQ | ✅ 已下载到 NAS（16 文件，8 样本 × paired-end） |
| 参考基因组 | ✅ 已拷贝到 `data/reference/epicauta/` |
| `data/fastq/epicauta/` | ✅ 目录已创建 |
| samples.tsv | ❌ 未创建 |
| config_epicauta.yaml | ❌ 未创建 |

### 样本映射（通过 spot 计数精确匹配）

| NAS 文件名 | Spots | SRA 描述 | 组 | 重命名为 | 
|-----------|:-----:|------|:--:|------|
| SRP293628.1 | 25.8M | DT rep2 | 🔴 跳过 | — | 
| SRP293628.2 | 23.9M | DT rep1 | 🔴 跳过 | — |
| SRP293628.3 | 26.1M | D rep3 | diapause | D3 |
| SRP293628.4 | 21.0M | D rep2 | diapause | D2 |
| SRP293628.5 | 21.9M | D rep1 | diapause | D1 |
| SRP293628.6 | 20.7M | ND rep3 | non-diapause | ND3 |
| SRP293628.7 | 24.1M | ND rep2 | non-diapause | ND2 |
| SRP293628.8 | 28.6M | ND rep1 | non-diapause | ND1 |

SRA 确认: SRP293628 包含 9 个样本（ND 3 + D 3 + DT 3），NAS 上下载了 8 个（缺 DT rep3, SRR13102881）。

### 实验设计

| 组 | 样本 | NAS 源文件 | 重复数 | 平均 spots |
|----|------|-----------|:-----:|:----------:|
| non-diapause | ND1, ND2, ND3 | SRP293628.8/.7/.6 | 3 | 24.5M |
| diapause | D1, D2, D3 | SRP293628.5/.4/.3 | 3 | 23.0M |

对比: `diapause_vs_non-diapause`  


文库: HiSeq X Ten, paired-end, 150bp, RNA-seq (strategy: other)

---

## Work Objectives

### Core Objective
在 Epicauta impressicornis 上运行 OmniQuant-RNA 两遍（decontam ON/OFF），产出共识 DEA + decontam 对比数据。

### Concrete Deliverables
- `data/fastq/epicauta/` 下 12 个 FASTQ 符号链接（6 样本 × paired-end）
- `data/fastq/epicauta/samples.tsv`
- `config/config_epicauta.yaml`
- `results.backup.epicauta_on/` — decontam ON 完整结果
- `results.backup.epicauta_off/` — decontam OFF 完整结果
- `results/benchmark/epicauta/` — decontam COMPARE 输出
- 论文 §3.6 填入核心数字

---

## Verification Strategy

- `ls data/fastq/epicauta/*_R1.fastq.gz | wc -l` = 6
- `snakemake --dry-run` 无错误
- `results.backup.epicauta_on/07.consensus_expression/diapause_vs_non-diapause/consensus_results.tsv` 存在且 >100 行
- OFF 同理
- `results/benchmark/epicauta/decontam_comparison_summary.tsv` 存在

---

## TODOs

- [x] 1. 创建 FASTQ 符号链接到 `data/fastq/epicauta/`

  **What to do**:
  - 从 NAS 创建符号链接，按映射表重命名：
    ```
    SRP293628.8_1.fastq.gz → ND1_R1.fastq.gz    SRP293628.8_2.fastq.gz → ND1_R2.fastq.gz
    SRP293628.7_1.fastq.gz → ND2_R1.fastq.gz    SRP293628.7_2.fastq.gz → ND2_R2.fastq.gz
    SRP293628.6_1.fastq.gz → ND3_R1.fastq.gz    SRP293628.6_2.fastq.gz → ND3_R2.fastq.gz
    SRP293628.5_1.fastq.gz → D1_R1.fastq.gz     SRP293628.5_2.fastq.gz → D1_R2.fastq.gz
    SRP293628.4_1.fastq.gz → D2_R1.fastq.gz     SRP293628.4_2.fastq.gz → D2_R2.fastq.gz
    SRP293628.3_1.fastq.gz → D3_R1.fastq.gz     SRP293628.3_2.fastq.gz → D3_R2.fastq.gz
    ```
  - 跳过 SRP293628.1 和 SRP293628.2（DT 组）
  - 使用符号链接（`ln -s`）而非复制，节省磁盘空间

  **Must NOT do**:
  - 不要复制文件（NAS 上已有 16 个 × ~2GB = ~32GB）
  - 不要包含 DT 组

  **QA Scenarios**:
  ```
  Scenario: 12 个符号链接正确创建
    Tool: bash
    Steps: ls -la data/fastq/epicauta/ND*_R1.fastq.gz data/fastq/epicauta/D*_R1.fastq.gz
    Expected: 6 files, all symlinks pointing to NAS paths
    Evidence: .sisyphus/evidence/epicauta-symlinks.txt
  ```

  **Commit**: NO

- [x] 2. 创建 `samples.tsv` 并备份当前 config

  **What to do**:
  - 创建 `data/fastq/epicauta/samples.tsv`:
    ```
    sample	fq1	fq2	group
    ND1	data/fastq/epicauta/ND1_R1.fastq.gz	data/fastq/epicauta/ND1_R2.fastq.gz	non-diapause
    ND2	data/fastq/epicauta/ND2_R1.fastq.gz	data/fastq/epicauta/ND2_R2.fastq.gz	non-diapause
    ND3	data/fastq/epicauta/ND3_R1.fastq.gz	data/fastq/epicauta/ND3_R2.fastq.gz	non-diapause
    D1	data/fastq/epicauta/D1_R1.fastq.gz	data/fastq/epicauta/D1_R2.fastq.gz	diapause
    D2	data/fastq/epicauta/D2_R1.fastq.gz	data/fastq/epicauta/D2_R2.fastq.gz	diapause
    D3	data/fastq/epicauta/D3_R1.fastq.gz	data/fastq/epicauta/D3_R2.fastq.gz	diapause
    ```
  - 备份当前 config: `cp config/config.yaml config/config_drosophila.yaml`
  - 确认参考基因组: `ls -la data/reference/epicauta/genome.fasta data/reference/epicauta/annotation.gff3`

  **QA Scenarios**:
  ```
  Scenario: samples.tsv 格式正确
    Tool: bash
    Steps: head -1 data/fastq/epicauta/samples.tsv && wc -l data/fastq/epicauta/samples.tsv
    Expected: header "sample	fq1	fq2	group", 7 lines (1 header + 6 samples)
    Evidence: .sisyphus/evidence/epicauta-samples-check.txt
  ```

  **Commit**: NO

- [x] 3. 创建 Epicauta 专用 config

  **What to do**:
  - 基于 `config/config.yaml` 创建 `config/config_epicauta.yaml`
  - 修改以下项:
    ```yaml
    samples: "data/fastq/epicauta/samples.tsv"
    samples_path: "data/fastq/epicauta"
    read_length: 150

    reference:
      genome: "data/reference/epicauta/genome.fasta"
      annotation: "data/reference/epicauta/annotation.gff3"

    decontam:
      enabled: true
      references:
        host_genome: "data/reference/epicauta/genome.fasta"
        host_transcriptome: "data/reference/epicauta/transcriptome.fasta"  # 自动生成
      # classifier db, ercc, technical_contam 等保持原路径

    dea:
      comparisons:
        - "diapause_vs_non-diapause"
    ```
  - 验证 NAS 路径可访问: `ls /mnt/nas/Database/k2_pluspfp_08_GB_20260226/`
  - 其他参数保持不变（fastp, hisat2, 定量器, consensus tiers 等）

  **QA Scenarios**:
  ```
  Scenario: config YAML 语法正确
    Tool: bash
    Steps: python3 -c "import yaml; c=yaml.safe_load(open('config/config_epicauta.yaml')); print('OK:', c['samples'])"
    Expected: "OK: data/fastq/epicauta/samples.tsv"
    Evidence: .sisyphus/evidence/epicauta-config-check.txt
  ```

  **Commit**: NO

- [ ] 4. 跑 decontam ON 管道 + 备份结果

  **What to do**:
  - 确认 `results/` 不存在或为空: `ls results/ 2>/dev/null || echo "results/ 不存在（OK）"`
  - 如果有残留: `rm -rf results/`
  - 切换 config: `cp config/config_epicauta.yaml config/config.yaml`
  - 保存运行时 config: `cp config/config.yaml config/config_epicauta_on.yaml`
  - 确认 `decontam.enabled: true`
  - 干运行检查: `snakemake --dry-run 2>&1 | tail -5`
  - 正式运行: `snakemake --use-conda --cores 32 --rerun-incomplete 2>&1 | tee logs/epicauta_on.log`
  - **备份结果**:
    ```bash
    mv results results.backup.epicauta_on
    echo "备份完成，验证目录:" && ls results.backup.epicauta_on/
    ```
  - 验证关键输出存在:
    - `results.backup.epicauta_on/03.decontam/stats/project_decontam_summary.tsv`
    - `results.backup.epicauta_on/07.consensus_expression/diapause_vs_non-diapause/consensus_results.tsv`
    - `results.backup.epicauta_on/03.decontam/clues/tables/sample_microbial_burden.tsv`
    - `results.backup.epicauta_on/08.reports/multiqc_report.html`

  **Must NOT do**:
  - 不要在 `results/` 已有内容时直接跑
  - 不要中断运行
  - 备份完成后不要删除 `results.backup.epicauta_on/`

  **QA Scenarios**:
  ```
  Scenario: ON 管道完成且备份包含全部 8 个 stage
    Tool: bash
    Steps: ls results.backup.epicauta_on/ | wc -l
    Expected: >= 8 (00-08 stages + 可能 benchmark)
    Evidence: .sisyphus/evidence/epicauta-on-complete.txt

  Scenario: consensus_results.tsv 存在且非空
    Tool: bash
    Steps: wc -l results.backup.epicauta_on/07.consensus_expression/diapause_vs_non-diapause/consensus_results.tsv
    Expected: > 100 行
    Evidence: .sisyphus/evidence/epicauta-on-consensus.txt
  ```

  **Commit**: NO

- [ ] 5. 跑 decontam OFF 管道 + 备份结果

  **What to do**:
  - 确认 `results/` 不存在或为空: `ls results/ 2>/dev/null || echo "results/ 不存在（OK）"`
  - 从 ON config 改 decontam OFF:
    ```bash
    cp config/config_epicauta_on.yaml config/config.yaml
    # 编辑 config.yaml: decontam.enabled: false
    ```
  - 保存运行时 config: `cp config/config.yaml config/config_epicauta_off.yaml`
  - 正式运行: `snakemake --use-conda --cores 32 --rerun-incomplete 2>&1 | tee logs/epicauta_off.log`
  - **备份结果**:
    ```bash
    mv results results.backup.epicauta_off
    echo "备份完成，验证目录:" && ls results.backup.epicauta_off/
    ```
  - 验证关键输出存在:
    - `results.backup.epicauta_off/07.consensus_expression/diapause_vs_non-diapause/consensus_results.tsv`
    - `results.backup.epicauta_off/08.reports/multiqc_report.html`
  - 恢复 config 到 ON 状态: `cp config/config_epicauta_on.yaml config/config.yaml`

  **Must NOT do**:
  - 不要在 `results/` 已有内容时直接跑
  - 不要用同一份 results/ 跑两次（ON 和 OFF 必须先后分开跑）

  **QA Scenarios**:
  ```
  Scenario: OFF 管道完成且备份包含 consensus
    Tool: bash
    Steps: ls results.backup.epicauta_off/ | wc -l
    Expected: >= 8 stages
    Evidence: .sisyphus/evidence/epicauta-off-complete.txt

  Scenario: ON 和 OFF 的 consensus 文件都独立存在
    Tool: bash
    Steps: ls results.backup.epicauta_on/07.consensus_expression/diapause_vs_non-diapause/consensus_results.tsv && ls results.backup.epicauta_off/07.consensus_expression/diapause_vs_non-diapause/consensus_results.tsv
    Expected: both files exist
    Evidence: .sisyphus/evidence/epicauta-both-backups.txt
  ```

  **QA Scenarios**:
  ```
  Scenario: OFF 管道完成
    Tool: bash
    Steps: wc -l results.backup.epicauta_off/07.consensus_expression/diapause_vs_non-diapause/consensus_results.tsv
    Expected: > 100 行
    Evidence: .sisyphus/evidence/epicauta-off-complete.txt
  ```

  **Commit**: NO

- [ ] 6. 运行 decontam COMPARE

  **What to do**:
  ```bash
  Rscript workflow/scripts/benchmark_decontam_spikein.R --mode compare \
    --on-dir results.backup.epicauta_on \
    --off-dir results.backup.epicauta_off \
    --contrast diapause_vs_non-diapause \
    --on-stats-dir results.backup.epicauta_on/03.decontam/stats \
    --off-stats-dir results.backup.epicauta_off/03.decontam/stats \
    --output-dir results/benchmark/epicauta
  ```
  - 检查输出: `ls results/benchmark/epicauta/decontam_comparison_summary.tsv`

  **QA Scenarios**:
  ```
  Scenario: COMPARE 输出包含关键指标
    Tool: bash
    Steps: grep "logFC_pearson_r\|tier_a_retention\|OFF_only" results/benchmark/epicauta/decontam_comparison_summary.tsv
    Expected: 3 lines with values
    Evidence: .sisyphus/evidence/epicauta-compare-complete.txt
  ```

  **Commit**: NO

- [ ] 7. 提取关键数字，更新论文 §3.6

  **What to do**:
  - 从 `consensus_summary.tsv` (ON) 提取:
    - universe_n, tier_a_n, tier_b_n, tier_c_n, conflict_n
  - 从 `decontam_comparison_summary.tsv` 提取:
    - logFC_pearson_r, tier_a_retention_pct, OFF_only_significant, ON_only_significant
  - 从 `sample_microbial_burden.tsv` 提取微生物负担概览
  - 从 `priority_targets.tsv` 提取重点目标（Wolbachia/病毒/真菌）存在性
  - 更新 `experiments/drafts/omniquant-paper-v1.md` §3.6:
    - Tier A 基因数 / 总基因宇宙
    - ON vs OFF decontam 关键数字（r, retention, OFF-only）
    - 微生物线索侧支的主要发现
    - 1-2 句生物学意义（滞育相关通路期望）

  **Must NOT do**:
  - 不要说 "ground truth" 或 "validation"
  - 不要和 Drosophila benchmark 比数字大小
  - 描述用 "demonstration of applicability" 语气

  **QA Scenarios**:
  ```
  Scenario: 论文 §3.6 包含具体数字
    Tool: bash
    Steps: grep -c "Tier A" experiments/drafts/omniquant-paper-v1.md
    Expected: >= 5 (至少 Abstract + Results + Epicauta 各提到)
    Evidence: .sisyphus/evidence/epicauta-paper-updated.txt
  ```

  **Commit**: NO

---

## Execution Strategy

```
Task 1 (符号链接) ── Task 2 (samples.tsv) ── Task 3 (config)
                                                  │
                        Task 4 (decontam ON) ◄────┘
                              │
                        Task 5 (decontam OFF)
                              │
                        Task 6 (COMPARE)
                              │
                        Task 7 (论文 §3.6)
```

串行执行，7 步。Task 1-3 为准备（~5 分钟），Task 4-5 为管道运行（~16 小时），Task 6-7 为分析（~15 分钟）。

## Commit Strategy

无 commit——所有产出均为 data/results/ 下的文件，按仓库规范不提交。

---

## Success Criteria

- [ ] 12 个符号链接正确指向 NAS FASTQ
- [ ] `snakemake --dry-run` 无错误
- [ ] decontam ON 管道输出包含 consensus_results.tsv（≥100 行）
- [ ] decontam OFF 管道输出包含 consensus_results.tsv（≥100 行）
- [ ] decontam COMPARE 输出包含 logFC_pearson_r 和 tier_a_retention
- [ ] 论文 §3.6 填入 Tier A 基因数 + 去污染关键数字 + 微生物负担概览
