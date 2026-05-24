# Reference 目录重构方案

## TL;DR

> **目标**: `data/reference/` 从"多基因组扁平混放 + 扁平派生文件互相覆盖"改为 `data/reference/{species}/` 隔离结构，每个物种自包含输入 + 派生文件 + 索引。
>
> **原则**: Config 一键切换 (`reference.species: drosophila_melanogaster`)，索引缓存隔离，不改动 Snakemake DAG 结构。
>
> **改动面**: 1 个核心 `.smk` 重写 + 7 个 `.smk` 机械替换 + 2 个 config yaml + 数据迁移。约 30 处硬编码字符串 → 变量引用。

---

## 1. 现状 vs 目标

### 当前（混乱）

```
data/reference/
├── Drosophila_melanogaster.BDGP6.54.dna.toplevel.fa    ← Drosophila 扁平
├── Drosophila_melanogaster.BDGP6.54.115.chr.gff3       ← Drosophila 扁平
├── epicauta/
│   ├── genome.fasta                                     ← Epicauta 子目录
│   └── annotation.gff3                                  ← Epicauta 子目录
├── transcriptome.fasta                          ← 派生，切换物种覆盖！
├── genome.gtf                                   ← 派生，切换物种覆盖！
├── genome.featurecounts.gtf                     ← 派生，切换物种覆盖！
├── genome.featurecounts.summary.tsv
├── hisat2_index/                                ← 索引，切换物种覆盖！
├── salmon_index/                                ← 774MB，切换物种覆盖！
├── kallisto_index/                              ← 索引，切换物种覆盖！
└── annotation_conversion_complete.flag          ← Flag，切换物种覆盖！
```

**Config 引用碎片化**：
- Drosophila ON config: `reference.genome: data/reference/Drosophila_melanogaster.BDGP6.54.dna.toplevel.fa`
- Epicauta ON config: `reference.genome: data/reference/epicauta/genome.fasta`
- 两个 config 的命名规则完全不同，无法自动化

### 目标（隔离）

```
data/reference/
├── drosophila_melanogaster/
│   ├── genome.fasta                              ← symlink → 原始文件（或直接复制）
│   ├── annotation.gff3                          ← symlink → 原始文件
│   ├── transcriptome.fasta                      ← 派生，不互相覆盖
│   ├── genome.gtf                               ← 派生，不互相覆盖
│   ├── genome.featurecounts.gtf                 ← 派生
│   ├── genome.featurecounts.summary.tsv
│   ├── annotation_conversion_complete.flag
│   ├── hisat2_index/                            ← 301MB，缓存
│   ├── salmon_index/                            ← 774MB，缓存
│   └── kallisto_index/                          ← 15MB，缓存
│
├── epicauta_diapause/                           ← 重命名自 epicauta/
│   ├── genome.fasta                              ← 已有
│   ├── annotation.gff3                           ← 已有
│   └── ...（同上，跑管道时生成）
│
└── bombyx_mori/                                 ← 未来
    └── （待下载）
```

**Config 统一为一个键**：

```yaml
# 之前
reference:
  genome: data/reference/Drosophila_melanogaster.BDGP6.54.dna.toplevel.fa
  annotation: data/reference/Drosophila_melanogaster.BDGP6.54.115.chr.gff3

# 之后
reference:
  species: drosophila_melanogaster   # 唯一需要改的键
```

---

## 2. 核心架构改动：`reference_config.smk`

### 2.1 当前逻辑（问题所在）

```python
REFERENCE_DIR = "data/reference"                                    # 硬编码

DEFAULT_REFERENCE_GENOME = "data/reference/genome.fasta"           # 假设扁平
DEFAULT_REFERENCE_ANNOTATION = "data/reference/genome.gff3"        # 假设扁平
DEFAULT_REFERENCE_GTF = "data/reference/genome.gtf"                # 硬编码
DEFAULT_REFERENCE_GFF3 = "data/reference/genome.gff3"              # 硬编码
DEFAULT_REFERENCE_TRANSCRIPTOME = "data/reference/transcriptome.fasta"  # 硬编码

def auto_detect_references():
    files = os.listdir(REFERENCE_DIR)    # 只扫描顶层，不递归！
    # ... 找到第一个 .fa 和第一个 .gff3 就用
```

### 2.2 新逻辑

```python
# === 从 config 读取物种名，派生所有路径 ===
REFERENCE_SPECIES = config["reference"]["species"]
REFERENCE_DIR = f"data/reference/{REFERENCE_SPECIES}"

REFERENCE_GENOME = f"{REFERENCE_DIR}/genome.fasta"
REFERENCE_ANNOTATION = f"{REFERENCE_DIR}/annotation.gff3"
REFERENCE_GTF = f"{REFERENCE_DIR}/genome.gtf"
FEATURECOUNTS_GTF = f"{REFERENCE_DIR}/genome.featurecounts.gtf"
TRANSCRIPTOME_FASTA = f"{REFERENCE_DIR}/transcriptome.fasta"

HISAT2_INDEX_DIR = f"{REFERENCE_DIR}/hisat2_index"
SALMON_INDEX_DIR = f"{REFERENCE_DIR}/salmon_index"
KALLISTO_INDEX_DIR = f"{REFERENCE_DIR}/kallisto_index"

# === 严格约定检查（解析阶段 fail fast） ===
if not os.path.isfile(REFERENCE_GENOME):
    raise ValueError(f"Missing genome fasta: {REFERENCE_GENOME}")
if not os.path.isfile(REFERENCE_ANNOTATION):
    raise ValueError(f"Missing annotation gff3: {REFERENCE_ANNOTATION}")

# annotation 格式检测
REFERENCE_SOURCE_ANNOTATION = REFERENCE_ANNOTATION
REFERENCE_TRANSCRIPTOME = TRANSCRIPTOME_FASTA

# 如果 annotation 本身就是 GTF（非 gff3），直接用，无需转换
if REFERENCE_ANNOTATION.endswith(".gtf"):
    REFERENCE_GTF = REFERENCE_ANNOTATION
    REFERENCE_SOURCE_FORMAT = "gtf"
else:
    REFERENCE_SOURCE_FORMAT = "gff3"
```

### 2.3 变化对比

| 变量 | 旧值 | 新值 |
|------|------|------|
| `REFERENCE_DIR` | `"data/reference"` | `f"data/reference/{species}"` |
| `REFERENCE_GENOME` | config 或 `genome.fasta` | 固定 `{REFERENCE_DIR}/genome.fasta` |
| `REFERENCE_ANNOTATION` | config 或自动检测 | 固定 `{REFERENCE_DIR}/annotation.gff3` |
| `REFERENCE_TRANSCRIPTOME` | 硬编码扁平 | `{REFERENCE_DIR}/transcriptome.fasta` |
| `HISAT2_INDEX_DIR` | 不存在 | 新增 `{REFERENCE_DIR}/hisat2_index` |
| `SALMON_INDEX_DIR` | 不存在 | 新增 `{REFERENCE_DIR}/salmon_index` |
| `KALLISTO_INDEX_DIR` | 不存在 | 新增 `{REFERENCE_DIR}/kallisto_index` |
| `FEATURECOUNTS_GTF` | 硬编码扁平 | `{REFERENCE_DIR}/genome.featurecounts.gtf` |
| `TRANSCRIPTOME_FASTA` | 不存在 | `{REFERENCE_DIR}/transcriptome.fasta` |

**不再需要自动检测**。每个物种目录必须严格遵循：
```
data/reference/{species}/
├── genome.fasta          ← 必须存在
├── annotation.gff3       ← 必须存在
├── genome.gtf            ← 管道生成
├── genome.featurecounts.gtf  ← 管道生成
├── transcriptome.fasta   ← 管道生成
├── hisat2_index/         ← 管道生成
├── salmon_index/         ← 管道生成
└── kallisto_index/       ← 管道生成
```
缺失 `genome.fasta` 或 `annotation.gff3` → 解析阶段 `ValueError` 终止。

### 2.4 Decontam 路径简化（零代码改动）

`decontam.smk` 的 `get_decontam_reference(key, default)` 已有 fallback：

```python
# decontam.smk:58
genome = get_decontam_reference("host_genome", REFERENCE_GENOME)
# → config 未设 decontam.references.host_genome → 自动用 REFERENCE_GENOME

# decontam.smk:113
transcriptome = get_decontam_reference("host_transcriptome", REFERENCE_TRANSCRIPTOME)
# → config 未设 decontam.references.host_transcriptome → 自动用 REFERENCE_TRANSCRIPTOME
```

`REFERENCE_GENOME` / `REFERENCE_TRANSCRIPTOME` 已从 `reference.species` 派生 → decontam 自动跟随。**Config 中删除冗余行即生效**。

---

## 3. SMK 文件改动清单（机械替换）

总览：约 **25 处硬编码字符串 → 变量引用**，分布在 7 个文件中。

### 3.1 `Snakefile`

| 行 | 旧 | 新 |
|----|----|----|
| 35 | `"data/reference/annotation_conversion_complete.flag"` | `f"{REFERENCE_DIR}/annotation_conversion_complete.flag"` |
| 38 | `REFERENCE_TRANSCRIPTOME` | 不变（已经是变量） |
| 56 | `"data/reference/kallisto_index/transcriptome.idx"` | `f"{KALLISTO_INDEX_DIR}/transcriptome.idx"` |
| 57 | `"data/reference/salmon_index"` | `SALMON_INDEX_DIR` |
| 58 | `"data/reference/hisat2_index"` | `HISAT2_INDEX_DIR` |
| 79 | `"data/reference/genome.featurecounts.gtf"` | `FEATURECOUNTS_GTF` |
| 80 | `"data/reference/genome.featurecounts.summary.tsv"` | `f"{REFERENCE_DIR}/genome.featurecounts.summary.tsv"` |
| 153 | `"data/reference/kallisto_index/transcriptome.idx"` | `f"{KALLISTO_INDEX_DIR}/transcriptome.idx"` |
| 161 | `"data/reference/salmon_index"` | `SALMON_INDEX_DIR` |
| 176 | `"data/reference/genome.featurecounts.gtf"` | `FEATURECOUNTS_GTF` |
| 177 | `"data/reference/genome.featurecounts.summary.tsv"` | `f"{REFERENCE_DIR}/genome.featurecounts.summary.tsv"` |

### 3.2 `workflow/rules/alignment.smk`

| 行 | 旧 | 新 |
|----|----|----|
| 15 | `directory("data/reference/hisat2_index")` | `directory(HISAT2_INDEX_DIR)` |
| 37 | `index="data/reference/hisat2_index"` | `index=HISAT2_INDEX_DIR` |
| 65 | `directory("data/reference/star_index")` | `directory(f"{REFERENCE_DIR}/star_index")` |
| 94 | `index="data/reference/star_index"` | `index=f"{REFERENCE_DIR}/star_index"` |

### 3.3 `workflow/rules/annotation_conversion.smk`

| 行 | 旧 | 新 |
|----|----|----|
| 9 | `"data/reference/annotation_conversion_complete.flag"` | `f"{REFERENCE_DIR}/annotation_conversion_complete.flag"` |
| 32 | `"data/reference/annotation_conversion_complete.flag"` | `f"{REFERENCE_DIR}/annotation_conversion_complete.flag"` |

### 3.4 `workflow/rules/quantification_kallisto.smk`

| 行 | 旧 | 新 |
|----|----|----|
| 14 | `"data/reference/kallisto_index/transcriptome.idx"` | `f"{KALLISTO_INDEX_DIR}/transcriptome.idx"` |
| 27 | `"data/reference/kallisto_index/transcriptome.idx"` | `f"{KALLISTO_INDEX_DIR}/transcriptome.idx"` |

### 3.5 `workflow/rules/quantification_salmon.smk`

| 行 | 旧 | 新 |
|----|----|----|
| 15 | `directory("data/reference/salmon_index")` | `directory(SALMON_INDEX_DIR)` |
| 37 | `"data/reference/salmon_index"` | `SALMON_INDEX_DIR` |

### 3.6 `workflow/rules/quantification_featurecounts.smk`

| 行 | 旧 | 新 |
|----|----|----|
| 4 | `"data/reference/genome.featurecounts.gtf"` (fallback) | `FEATURECOUNTS_GTF`（变量已从 config 派生，fallback 用 `FEATURECOUNTS_GTF`） |
| 7 | `"data/reference/genome.featurecounts.summary.tsv"` (fallback) | `f"{REFERENCE_DIR}/genome.featurecounts.summary.tsv"` |

### 3.7 `workflow/modular_workflows/annotation_conversion_only.smk`

| 行 | 旧 | 新 |
|----|----|----|
| 30 | `"data/reference/annotation_conversion_complete.flag"` | `f"{REFERENCE_DIR}/annotation_conversion_complete.flag"` |
| 35 | `"data/reference/annotation_conversion_complete.flag"` | `f"{REFERENCE_DIR}/annotation_conversion_complete.flag"` |

### 3.8 `workflow/rules/benchmark.smk`

| 行 | 旧 | 新 |
|----|----|----|
| 7 | `BENCH_REFERENCE_DIR = "data/reference"` | `BENCH_REFERENCE_DIR = REFERENCE_DIR` |

### 3.9 不需要改动的文件

- **`decontam.smk`**：用 `get_decontam_reference()` 从 config 读路径，不硬编码 `data/reference/` ✅
- **`reference_namespace.smk`**：通过 `REFERENCE_GTF`、`REFERENCE_TX2GENE` 等变量引用，ok ✅
- **所有 Python 脚本**：通过 `snakemake.params` 接收路径，不硬编码 ✅
- **`common.smk`**：只 include `reference_config.smk`，无引用路径 ✅

---

## 4. Config 改动

### 4.1 `config/runs/drosophila_wolbachia/2026-05-24_decontam-kaiju.yaml`

```yaml
# 删除
reference:
  annotation: data/reference/Drosophila_melanogaster.BDGP6.54.115.chr.gff3
  genome: data/reference/Drosophila_melanogaster.BDGP6.54.dna.toplevel.fa

# 替换为
reference:
  species: drosophila_melanogaster
```

**decontam 引用无需改动**：`get_decontam_reference("host_genome", REFERENCE_GENOME)` 的 fallback 已经是 `REFERENCE_GENOME`（即 `data/reference/{species}/genome.fasta`）。如果 config 中 `decontam.references.host_genome` 未设置，自动使用 REFERENCE_GENOME。同理 `host_transcriptome` 自动使用 `REFERENCE_TRANSCRIPTOME`。

→ **删除** config 中冗余的 `decontam.references.host_genome` 和 `host_transcriptome` 行（如果存在），让 fallback 生效。

### 4.2 `config/runs/epicauta_diapause/2026-05-24_decontam-kaiju.yaml`

```yaml
# 删除
reference:
  annotation: data/reference/epicauta/annotation.gff3
  genome: data/reference/epicauta/genome.fasta

# 替换为
reference:
  species: epicauta_diapause
```

同上，删除 `decontam.references.host_genome`/`host_transcriptome` 行。

### 4.3 `config/config.yaml`（软链接）

当前指向 `runs/epicauta_diapause/2026-05-24_decontam-kaiju.yaml` → 自动生效。

---

## 5. 数据迁移步骤

### Step A: 创建物种子目录 + symlink 原始文件

```bash
# Drosophila — 使用绝对路径确保容器/集群兼容
mkdir -p data/reference/drosophila_melanogaster
ln -sf "$(realpath data/reference/Drosophila_melanogaster.BDGP6.54.dna.toplevel.fa)" \
       data/reference/drosophila_melanogaster/genome.fasta
ln -sf "$(realpath data/reference/Drosophila_melanogaster.BDGP6.54.115.chr.gff3)" \
       data/reference/drosophila_melanogaster/annotation.gff3

# Epicauta（重命名目录）
mv data/reference/epicauta data/reference/epicauta_diapause

# Bombyx（预留）
mkdir -p data/reference/bombyx_mori
```

> ⚠️ **Symlink 注意事项**：使用 `realpath` 解析为绝对路径，避免相对路径 `../` 在 Singularity/Apptainer 容器或集群环境中挂载失败。HISAT2/Salmon/kallisto 均支持读取 symlink 指向的实际文件。

### Step B: 删除旧扁平派生文件

```bash
rm -f data/reference/transcriptome.fasta
rm -f data/reference/genome.gtf
rm -f data/reference/genome.featurecounts.gtf
rm -f data/reference/genome.featurecounts.summary.tsv
rm -f data/reference/annotation_conversion_complete.flag
rm -rf data/reference/hisat2_index/
rm -rf data/reference/salmon_index/
rm -rf data/reference/kallisto_index/
```

### Step C: 重新运行管道生成新索引

切换 species 后首次运行会自动生成该物种的派生文件 + 索引。之后切换回其他物种不会覆盖。

```bash
# 当前 config.yaml → epicauta_diapause
snakemake --use-conda --cores 32   # 生成 epicauta 的所有派生文件 + 索引

# 切换到 Drosophila
ln -sf runs/drosophila_wolbachia/2026-05-24_decontam-kaiju.yaml config/config.yaml
snakemake --use-conda --cores 32   # 生成 drosophila 的所有派生文件 + 索引
# ← epicauta 的索引不会被覆盖！
```

---

## 6. 验证清单

跑完代码改动后：

- [ ] `snakemake --dry-run` 无报错，DAG 完整
- [ ] `grep -rn 'data/reference/[^{]' workflow/rules/ Snakefile workflow/modular_workflows/` 无残留硬编码
- [ ] 切换 `reference.species` → 清除 `results/` → 重新跑管道 → 输出正确
- [ ] 两个物种的索引目录共存，`du -sh data/reference/*/hisat2_index` 显示两个独立索引
- [ ] Config 中 `decontam.references.host_genome` 已删除，fallback 生效
- [ ] `config/runs/` 下的两个 yaml 只有 `reference.species` 一个键
- [ ] 缺失 `genome.fasta` 或 `annotation.gff3` 时，`snakemake --dry-run` 直接 ValueError 终止（fail fast）

---

## 7. 风险 & 回滚

| 风险 | 缓解 |
|------|------|
| 原始 Drosophila 文件保留在根目录 | 用 symlink 而非移动，原始文件不动 |
| Symlink 在容器/集群中失效 | 使用 `realpath` 解析为**绝对路径**。HISAT2/Salmon/kallisto 均支持读取 symlink |
| `bombyx_mori` 空目录导致报错 | 预期行为：fail fast。放好 genome.fasta + annotation.gff3 后自动通过 |
| decontam 路径重复配置 | 已精简：`get_decontam_reference()` fallback 自动使用 `REFERENCE_GENOME`，config 中删除冗余行 |
| benchmark.smk 用了 `REFERENCE_DIR` | benchmark.smk 是调试用，不影响主线 |
| 改动面大（7 个 SMK 文件） | 机械替换，`grep` + `snakemake --dry-run` 一次性验证 |
| `reference.genome` / `reference.annotation` 被删后旧脚本引用 | 不保留向后兼容。`reference_config.smk` 直接从 `reference.species` 派生所有路径 |

---

*Plan v1. 2026-05-24. 基于 explore agents 对 25 处硬编码路径的全面审计。*
