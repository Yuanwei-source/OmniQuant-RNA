# OmniQuant-RNA

OmniQuant-RNA 是一个基于 Snakemake 的 RNA-seq 自动化分析流程，面向多定量器并行分析、差异表达统计整合，以及在统一基因命名空间下进行共识差异表达分析。当前主流程已经覆盖原始数据质控、修剪、比对、定量、差异表达和共识输出，并支持在非模式昆虫等注释不完备场景下进行较稳健的 gene-level 整合分析。

## 主要能力

- 多定量器并行：featureCounts、StringTie、Salmon、Kallisto。
- 差异表达分析：DESeq2、edgeR、limma。
- 共识差异表达：基于多定量器 DESeq2 结果进行 RRA 和 CCT 聚合。
- 统一命名空间：构建 reference tx2gene、gene namespace、StringTie bridge 和 master tx2gene。
- 参考自动识别：自动接管参考基因组和注释文件，并生成所需索引与派生注释。
- 模块化运行：支持完整流程和按模块执行。
- 日志与结果分层：定量原生输出、正式矩阵、审计表和共识诊断分开组织。

## 工作流概览

主流程按以下顺序执行：

1. 参考文件标准化与注释转换。
2. 原始 reads 质控和 fastp 修剪。
3. HISAT2 比对与 BAM 生成。
4. featureCounts、StringTie、Salmon、Kallisto 定量。
5. 各定量器独立差异表达分析。
6. 多定量器共识差异表达整合与可视化。
7. MultiQC 与汇总报告输出。

## 目录结构

```text
OmniQuant-RNA/
├── README.md
├── Snakefile
├── run_analysis.sh
├── run_modular.sh
├── config/
│   └── config.yaml
├── data/
│   ├── fastq/
│   │   └── samples.tsv
│   └── reference/
├── envs/
│   ├── alignment.yaml
│   ├── dea.yaml
│   ├── featurecounts.yaml
│   ├── qc.yaml
│   ├── quantification.yaml
│   └── stringtie.yaml
├── workflow/
│   ├── modular_workflows/
│   ├── rules/
│   └── scripts/
├── logs/
└── results/
    ├── 00.reference/
    ├── 01.raw_qc/
    ├── 02.trimmed_data/
    ├── 03.decontam/
    ├── 04.alignment/
    ├── 05.quantification/
    ├── 06.differential_expression/
    ├── 07.consensus_expression/
  └── 08.reports/
```

## 环境要求

- Linux 优先。
- Conda 或 Mamba。
- Snakemake。
- 足够的 CPU、内存和磁盘空间。

一个常见安装方式如下：

```bash
conda install -c conda-forge -c bioconda snakemake
```

## 输入准备

### 1. FASTQ 样本表

样本表位于 [data/fastq/samples.tsv](data/fastq/samples.tsv)，至少应包含以下列：

```text
sample        fq1                                fq2                                group
Control-1     /path/to/Control-1_R1.fastq.gz     /path/to/Control-1_R2.fastq.gz     Control
Control-2     /path/to/Control-2_R1.fastq.gz     /path/to/Control-2_R2.fastq.gz     Control
Treatment-1   /path/to/Treatment-1_R1.fastq.gz   /path/to/Treatment-1_R2.fastq.gz   Treatment
Treatment-2   /path/to/Treatment-2_R1.fastq.gz   /path/to/Treatment-2_R2.fastq.gz   Treatment
```

可使用以下脚本自动生成样本表：

```bash
python workflow/scripts/generate_samples.py ./data/fastq -o ./data/fastq/samples.tsv
```

当前 FastQC 输入兼容常见的 paired-end 命名模式，并会在内部生成标准化别名，确保输出稳定对应到 sample_R1 和 sample_R2。

### 2. 参考文件

将参考文件放入 [data/reference](data/reference)。流程会自动识别并标准化以下内容：

- 基因组 FASTA。
- GFF3 或 GTF 注释。
- 转录组 FASTA，若缺失可自动提取。

主配置见 [config/config.yaml](config/config.yaml)。

## 关键配置项

以下几组配置与当前主流程直接相关。

### 参考与样本

```yaml
samples: "data/fastq/samples.tsv"
read_length: 150

reference:
  genome: "data/reference/genome.fasta"
  gff3: "data/reference/genome.gff3"
  gtf: "data/reference/genome.gtf"
  transcriptome: "data/reference/transcriptome.fasta"
```

### 多定量器统一命名空间

```yaml
namespace:
  output_dir: "results/00.reference"
  tx2gene_reference: "results/00.reference/tx2gene_reference.tsv"
  stringtie_bridge: "results/00.reference/stringtie_tx2gene_bridge.tsv"
  tx2gene_master: "results/00.reference/tx2gene_master.tsv"
  gene_namespace: "results/00.reference/gene_namespace.tsv"
```

### DEA 导入策略

```yaml
dea_import:
  main_only: true
  mapping_policy: "conservative"
  novel_policy: "discovery_only"
```

### 共识差异表达

```yaml
consensus:
  enabled: true
  methods: ["deseq2"]
  quantifiers: ["featurecounts", "stringtie", "salmon", "kallisto"]
  contrasts: "all"
  p_clip: 1.0e-16
```

### 单定量器 DEA

```yaml
dea:
  methods: ["deseq2", "edger", "limma"]
  comparisons: "all"
  batch_column: null
  fdr_threshold: 0.05
  lfc_threshold: 1.0
  min_count: 10
  min_methods: 3
```

## 运行方式

### 完整流程

推荐使用入口脚本：

```bash
./run_analysis.sh
```

也可以直接运行 Snakemake：

```bash
snakemake --use-conda --cores 16 --rerun-incomplete
```

后台运行示例：

```bash
nohup ./run_analysis.sh > analysis.log 2>&1 &
```

### 干运行检查

```bash
snakemake --dry-run
snakemake --dry-run --quiet
snakemake --rulegraph | dot -Tpng > rules_graph.png
```

### 模块化运行

[run_modular.sh](run_modular.sh) 当前支持以下模块：

- `all`：运行完整工作流。
- `qc`：只运行质量控制。
- `quantification`：只运行定量阶段，输出到 [results/05.quantification](results/05.quantification)，并按 [results/05.quantification/native](results/05.quantification/native)、[results/05.quantification/matrices](results/05.quantification/matrices)、[results/05.quantification/audit](results/05.quantification/audit) 分层组织。
- `alignment`：只运行序列比对。
- `differential_expression` / `dea`：只运行差异表达分析，读取 [results/05.quantification/matrices](results/05.quantification/matrices) 和 [results/00.reference/import_manifests](results/00.reference/import_manifests) 作为正式输入。
- `dry-run`：干运行检查主流程及关键模块。

示例：

```bash
./run_modular.sh qc --cores 8
./run_modular.sh quantification --cores 16
./run_modular.sh dea --cores 16
./run_modular.sh all --cores 24
./run_modular.sh dry-run
```

对应的模块化 Snakefile 位于 [workflow/modular_workflows](workflow/modular_workflows)。

## 主要输出

### 1. 统一参考空间

位于 [results/00.reference](results/00.reference)：

- tx2gene_reference.tsv
- gene_namespace.tsv
- stringtie_tx2gene_bridge.tsv
- tx2gene_master.tsv
- import_manifests/

### 2. 质量控制

- [results/01.raw_qc](results/01.raw_qc)：原始 reads FastQC 结果。
- [results/02.trimmed_data](results/02.trimmed_data)：fastp 输出和修剪后 FastQC 结果。
- [results/03.decontam](results/03.decontam)：去污染 clean reads、统计表和 clean FastQC 结果。
- [results/08.reports](results/08.reports)：MultiQC 报告。

### 3. 比对与定量

- [results/04.alignment](results/04.alignment)：HISAT2 或 STAR 比对结果。
- [results/05.quantification](results/05.quantification)：定量阶段统一输出目录。
- 其中 [results/05.quantification/native](results/05.quantification/native) 保存各定量器原生结果。
- 其中 [results/05.quantification/matrices](results/05.quantification/matrices) 保存下游分析使用的正式矩阵。
- 其中 [results/05.quantification/audit](results/05.quantification/audit) 保存轻量审计表，如 StringTie gene_id 映射。

### 4. 差异表达分析

位于 [results/06.differential_expression](results/06.differential_expression)，每个 quantifier 目录下包含：

- 原始 DEA 结果表。
- dea_session.rds。
- normalized_counts.csv。
- integration 子目录中的整合图表。

### 5. 共识差异表达分析

位于 [results/07.consensus_expression](results/07.consensus_expression)，按 contrast 分目录输出：

- consensus_results.tsv
- consensus_summary.tsv
- tier_diagnostics.tsv
- significance_membership.tsv
- sensitivity_analysis.tsv
- logFC_scatter_salmon_vs_featurecounts.pdf
- consensus_volcano.pdf
- significance_upset.pdf

## 方法说明

### 单定量器层

- featureCounts 直接使用 gene-level count matrix。
- Salmon、Kallisto 和 StringTie 通过统一 tx2gene 体系导入 gene-level 分析空间。
- StringTie 采用保守映射策略，优先保留 reference-compatible 特征进入主分析。

### 共识层

共识分析当前以各定量器的 DESeq2 结果为输入，对共享基因宇宙进行整合，并输出分层证据结果。

核心指标包括：

- support_n
- sign_consistency_n
- consensus_direction
- consensus_logFC
- logFC_CV
- best_rra_fdr
- best_cct_fdr

默认将基因划分为 Tier A、Tier B、Tier C 和 unclassified。

## 日志与调试

- Snakemake 主日志位于 [.snakemake/log](.snakemake/log)。
- 各规则日志位于 [logs](logs)。
- StringTie 聚合和 gene mapping 的详细 verbose 输出会写入对应日志文件，而不是直接打印到屏幕。

## 当前默认假设

- 当前主流程默认使用 HISAT2。
- 共识层默认启用。
- 共识层默认整合四个 quantifier 的 DESeq2 输出。
- 若使用数值开头的分组标签，流程内部会自动使用安全名称构建设计矩阵，但保留原始比较标签用于输出。

## 建议使用方式

如果你的目标是稳定获得最终 gene-level 结论，建议直接跑完整流程或至少跑到 DEA 与 consensus 输出层，而不是只停留在单一定量器矩阵阶段。
