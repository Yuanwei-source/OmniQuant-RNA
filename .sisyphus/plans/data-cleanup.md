# 数据与结果整理方案 v2

## TL;DR

> **目标**: 按最终目录结构清理数据/config/results，建立可复现、可追溯的 benchmark 管理体系。
>
> **原则**: Snakemake 不改、config 不变、结果不覆盖、每条 run 自包含。

---

## 1. 最终目录结构

```
OmniQuant-RNA/
├── workflow/                    # CODE: Snakemake (不动)
├── envs/                        # CODE: conda yaml (不动)
│
├── config/                      # 配置
│   ├── config.yaml              #   软链接 → runs/.../xxx.yaml
│   └── runs/
│       ├── drosophila_wolbachia/
│       │   ├── 2026-05-24_decontam-kaiju.yaml
│       │   └── 2026-05-24_decontam-off.yaml
│       └── epicauta_diapause/
│           ├── 2026-05-24_decontam-kaiju.yaml
│           └── 2026-05-24_decontam-off.yaml
│
├── data/                        # 输入
│   ├── drosophila_wolbachia/
│   │   └── samples.tsv, *_R1.fastq.gz, *_R2.fastq.gz
│   ├── bombyx_mori/
│   └── epicauta_diapause/
│       └── samples.tsv, *_R1.fastq.gz→nas, *_R2.fastq.gz→nas
│
├── benchmark_results/           # 🆕 所有 benchmark 产出
│   ├── drosophila_wolbachia/
│   │   ├── runs/
│   │   │   ├── 2026-05-24_decontam-kaiju/
│   │   │   │   ├── 00.reference/
│   │   │   │   ├── ...
│   │   │   │   ├── 07.consensus_expression/
│   │   │   │   ├── RUN_METADATA.yaml
│   │   │   │   ├── config_snapshot.yaml
│   │   │   │   ├── samplesheet_snapshot.tsv
│   │   │   │   └── _RUN_COMPLETE
│   │   │   └── 2026-05-24_decontam-off/
│   │   ├── comparisons/
│   │   │   └── 2026-05-24_kaiju-vs-off/
│   │   │       ├── decontam_comparison_summary.tsv
│   │   │       └── figures/
│   │   └── archive/
│   │
│   ├── epicauta_diapause/
│   │   ├── runs/
│   │   ├── comparisons/
│   │   └── archive/
│   │
│   ├── bombyx_mori/
│   │   └── runs/
│   │
│   └── cross_species/                    # 跨数据集 benchmark
│       ├── annotation_degradation/
│       ├── subsampling_stability/
│       ├── ablation/
│       ├── figures/
│       ├── statistical_tests.tsv
│       └── FINAL_REPORT.md
│
├── results/                    # ⚠️ 临时 ! Snakemake 运行时存在，跑完即 mv
│
└── experiments/                # 论文
    └── drafts/
        └── omniquant-paper-v1.md
```

### 每条 run 自包含

```
benchmark_results/{dataset}/runs/{date}_{config}/
├── 00.reference/ ... 07.consensus_expression/    # 管道输出
├── RUN_METADATA.yaml          # 运行时间、数据集、对比、config 名
├── config_snapshot.yaml       # 跑时的 config 完整快照
├── samplesheet_snapshot.tsv   # 跑时的样本表
├── software_versions.tsv      # 工具版本（可选）
└── _RUN_COMPLETE              # 空文件，标记管道正常结束
```

---

## 2. 清理执行步骤

### Step 1: 清理过期/残留文件
```bash
# 下载残留
rm -f data/fastq/wolbachia/SRR14277097.sra.prf
rm -f data/fastq/wolbachia/SRR14277097.sra.tmp

# 过期备份
rm -f config/config.yaml.drosophila.bak

# 散落临时文件
rm -f data/fastq/samples.tsv          # 旧根目录样本表（已不用）
rm -f data/fastq/subsample.sh         # 一次性脚本
```

### Step 2: 重整 config 目录
```bash
mkdir -p config/runs/{drosophila_wolbachia,epicauta_diapause}

# 现有 ON/OFF config 移入 runs/
mv config/config_drosophila_on.yaml  config/runs/drosophila_wolbachia/2026-05-24_decontam-kaiju.yaml
mv config/config_epicauta_on.yaml    config/runs/epicauta_diapause/2026-05-24_decontam-kaiju.yaml

# 删掉多余副本
rm config/config_drosophila.yaml config/config_epicauta.yaml

# 设置软链接
ln -sf runs/epicauta_diapause/2026-05-24_decontam-kaiju.yaml config/config.yaml
```

### Step 3: 清理当前 results/
```bash
rm -rf results/ .snakemake/locks/
```

### Step 4: 整理旧 benchmark results
```bash
mkdir -p benchmark_results/{drosophila_wolbachia,epicauta_diapause,bombyx_mori}/archive
mkdir -p benchmark_results/{drosophila_wolbachia,epicauta_diapause}/{runs,comparisons}
mkdir -p benchmark_results/cross_species/{annotation_degradation,subsampling_stability,ablation,figures}

# 归档旧 Kraken2 结果
mv results.backup.decontam_on  benchmark_results/drosophila_wolbachia/archive/2026-05-08_decontam-kraken2/
mv results.backup.decontam_off benchmark_results/drosophila_wolbachia/archive/2026-05-08_decontam-off/

# 整理 Epicauta 旧结果
mv results.backup.epicauta_on_fullk2 benchmark_results/epicauta_diapause/archive/2026-05-20_decontam-kraken2-fullnt/
mv results.backup.epicauta_on_fixed   benchmark_results/epicauta_diapause/archive/2026-05-18_decontam-kraken2-pluspfp/
mv results.backup.epicauta_on         benchmark_results/epicauta_diapause/archive/2026-05-15_decontam-kraken2-pluspfp/
mv results.backup.epicauta_off        benchmark_results/epicauta_diapause/archive/2026-05-15_decontam-off/

# 跨物种 benchmark 产出
mv results/benchmark/*.tsv results/benchmark/figures/ benchmark_results/cross_species/ 2>/dev/null
```

### Step 5: 重整数据目录（⚠️ 手动确认后执行）
```bash
# 确认 sub_SRR14101759-64 归属后移入对应目录
# 确认 wolbachia/SRR14277094-99 是否需要删除并重新下载
```

---

## 3. 日常使用流程

```bash
# === 跑之前 ===
ln -sf runs/drosophila_wolbachia/2026-05-24_decontam-kaiju.yaml config/config.yaml
rm -rf results/ .snakemake/locks/

# === 跑 ===
snakemake --use-conda --cores 32 --rerun-incomplete

# === 跑完后 ===
RUN_DIR="benchmark_results/drosophila_wolbachia/runs/2026-05-24_decontam-kaiju"
mkdir -p "$RUN_DIR"
mv results/* "$RUN_DIR/"
touch "$RUN_DIR/_RUN_COMPLETE"
cp config/config.yaml "$RUN_DIR/config_snapshot.yaml"

# === COMPARE 分析 ===
Rscript workflow/scripts/benchmark_decontam_spikein.R --mode compare \
  --on-dir  benchmark_results/drosophila_wolbachia/runs/2026-05-24_decontam-kaiju \
  --off-dir benchmark_results/drosophila_wolbachia/runs/2026-05-24_decontam-off \
  --contrast 60d_vs_1d \
  --output-dir benchmark_results/drosophila_wolbachia/comparisons/2026-05-24_kaiju-vs-off
```

---

## 4. 待执行清单

| Step | 内容 | 破坏性 |
|:----:|------|:------:|
| 1 | 删除残留文件 | 否 |
| 2 | 重整 config 目录 | ⚠️ 路径变了 |
| 3 | 清理 results/ | 否 |
| 4 | 整理旧 benchmark results | ⚠️ 路径变了 |
| 5 | 重整数据目录 | ⚠️ 需确认归属 |
| — | 重新下载 Drosophila Wolbachia 数据 | ⚠️ SRA 下载 |

---

*Plan v2. 基于用户最终目录结构。*
