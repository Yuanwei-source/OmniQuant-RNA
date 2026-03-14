#!/bin/bash

if [ -z "${BASH_VERSION:-}" ]; then
    exec bash "$0" "$@"
fi

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
if [ ! -f "data/fastq/samples.tsv" ]; then
    echo "错误: 样本配置文件 data/fastq/samples.tsv 不存在"
    echo "请执行以下操作之一："
    echo "1. 手动创建该文件，格式需包含 sample 和 path 两列"
    echo "2. 运行脚本自动生成（需指定您的数据目录）："
    echo "   python workflow/scripts/generate_samples.py <您的fastq数据目录> -o data/fastq/samples.tsv"
    exit 1
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

# 询问用户是否继续 (如果在后台/nohup运行，则自动跳过询问并继续)
if [ -t 0 ]; then
    read -p "是否继续执行工作流? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "已取消执行"
        exit 0
    fi
else
    echo "检测到非交互环境(如nohup)，自动继续执行工作流..."
fi

# 运行工作流
echo "3. 开始运行工作流..."
CORES=${CORES:-32}  # 默认使用32个核心
snakemake --use-conda --cores $CORES --printshellcmds

echo "=== 工作流执行完成 ==="
echo "结束时间: $(date)"
echo ""
echo "主要结果文件:"
echo "- results/00.reference/: 统一 tx2gene / gene namespace / import manifest"
echo "- results/03.decontam/: 去污染 clean/stats/qc 输出"
echo "- results/05.quantification/native/: 各定量器原生输出"
echo "- results/05.quantification/matrices/: 正式矩阵输出"
echo "- results/05.quantification/audit/: 审计与映射辅助表"
echo "- results/06.differential_expression/: 正式 DEA 结果"
echo "- results/08.reports/multiqc_report.html: 综合质量控制报告"
