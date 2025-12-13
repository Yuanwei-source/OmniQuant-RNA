#!/bin/bash

# OmniQuant-RNA - 运行脚本
# 使用方法: ./run_analysis.sh

set -euo pipefail

# 检查 Snakemake 是否安装
if ! command -v snakemake &> /dev/null; then
    echo "错误: Snakemake 未安装. 请先安装 Snakemake:"
    echo "conda install -c bioconda snakemake"
    exit 1
fi

# 检查必要文件是否存在
if [ ! -f "config/samples.tsv" ]; then
    echo "警告: config/samples.tsv 文件不存在"
    echo "尝试自动生成样本文件..."
    if [ -d "data/fastq" ]; then
        python workflow/scripts/generate_samples.py data/fastq/ -o config/samples.tsv
        echo "已生成 config/samples.tsv"
    else
        echo "错误: 数据目录 data/fastq 不存在，无法生成样本文件"
        exit 1
    fi
fi

# 创建必要的目录
mkdir -p data/fastq data/reference

echo "=== OmniQuant-RNA 转录组定量分析工作流 ==="
echo "开始时间: $(date)"

# 检查工作流
echo "1. 检查工作流语法..."
snakemake --dry-run

# 显示将要执行的任务
echo "2. 将要执行的任务:"
snakemake --dry-run --quiet

# 询问用户是否继续
read -p "是否继续执行工作流? (y/N): " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "已取消执行"
    exit 0
fi

# 运行工作流
echo "3. 开始运行工作流..."
CORES=${CORES:-32}  # 默认使用32个核心
snakemake --use-conda --cores $CORES --printshellcmds

echo "=== 工作流执行完成 ==="
echo "结束时间: $(date)"
echo ""
echo "主要结果文件:"
echo "- results/multiqc_report.html: 综合质量控制报告"
echo "- results/quantification/: 各样本的定量结果"
echo "- results/fastqc/: 质量控制结果"
