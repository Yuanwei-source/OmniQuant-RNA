#!/bin/bash

# OmniQuant-RNA 模块化运行脚本
# 使用方法: ./run_modular_analysis.sh [module] [options]

set -euo pipefail

CORES=${CORES:-8}
CONDA_ARGS="--use-conda"

show_usage() {
    echo "OmniQuant-RNA 模块化分析工具"
    echo ""
    echo "使用方法: $0 [module] [options]"
    echo ""
    echo "可用模块:"
    echo "  all           - 运行完整工作流"
    echo "  qc            - 只运行质量控制"
    echo "  quantification- 只运行转录组定量"
    echo "  alignment     - 只运行序列比对"
    echo "  differential_expression - 只运行差异表达分析 (简写: dea)"
    echo "  dry-run       - 检查工作流（干运行）"
    echo ""
    echo "选项:"
    echo "  --cores N     - 使用 N 个CPU核心 (默认: 8)"
    echo "  --help        - 显示此帮助信息"
    echo ""
    echo "例子:"
    echo "  $0 qc                    # 只运行质量控制"
    echo "  $0 all --cores 16        # 使用16核运行完整流程"
    echo "  $0 dry-run              # 检查工作流"
}

# 解析命令行参数
MODULE=""
while [[ $# -gt 0 ]]; do
    case $1 in
        --cores)
            CORES="$2"
            shift 2
            ;;
        --help)
            show_usage
            exit 0
            ;;
        all|qc|quantification|alignment|differential_expression|dea|dry-run)
            MODULE="$1"
            shift
            ;;
        *)
            echo "未知参数: $1"
            show_usage
            exit 1
            ;;
    esac
done

if [[ -z "$MODULE" ]]; then
    echo "错误: 请指定要运行的模块"
    show_usage
    exit 1
fi

echo "=== OmniQuant-RNA 模块化分析 ==="
echo "模块: $MODULE"
echo "CPU核心数: $CORES"
echo "开始时间: $(date)"

case $MODULE in
    all)
        echo "运行完整工作流..."
        snakemake $CONDA_ARGS --cores $CORES --printshellcmds
        ;;
    qc)
        echo "运行质量控制模块..."
        snakemake -s workflow/qc_only.smk $CONDA_ARGS --cores $CORES --printshellcmds
        ;;
    quantification)
        echo "运行转录组定量模块..."
        snakemake -s workflow/quantification_only.smk $CONDA_ARGS --cores $CORES --printshellcmds
        ;;
    differential_expression|dea)
        echo "运行差异表达分析模块..."
        snakemake -s workflow/differential_expression_only.smk $CONDA_ARGS --cores $CORES --printshellcmds
        ;;
    alignment)
        echo "运行序列比对模块..."
        snakemake -s workflow/alignment_only.smk $CONDA_ARGS --cores $CORES --printshellcmds
        ;;
    dry-run)
        echo "检查工作流（干运行）..."
        snakemake --dry-run
        echo ""
        echo "各模块检查:"
        echo "- QC模块:"
        snakemake -s workflow/qc_only.smk --dry-run
        echo "- 定量模块:"
        snakemake -s workflow/quantification_only.smk --dry-run
        ;;
esac

echo "=== 分析完成 ==="
echo "结束时间: $(date)"
