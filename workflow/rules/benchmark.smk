# Benchmark rules — multi-quantifier consensus stability and annotation degradation
# These rules wrap standalone R benchmark scripts for one-command execution.

import os

# Local path constants — available in common.smk but referenced here for standalone use
BENCH_REFERENCE_DIR = REFERENCE_DIR

# ── Rule: subsampling stability benchmark ──────────────────────────────────────
rule benchmark_subsampling:
    """
    子采样稳定性基准：从 60d_vs_1d 对比中每组选 2 个重复，
    共 9 子集评估多定量器共识 DEA 的稳定性。
    """
    output:
        stability = "results/benchmark/subsampling_stability.tsv",
        summary   = "results/benchmark/subsampling_summary.tsv",
        fig_dist  = "results/benchmark/figures/stability_distribution.png",
        fig_curve = "results/benchmark/figures/stability_yield_curve.png"
    conda:
        "../../envs/dea.yaml"
    benchmark:
        "benchmarks/benchmark_subsampling.tsv"
    log:
        "logs/benchmark/subsampling.log"
    script:
        "../scripts/benchmark_subsampling.R"

# ── Rule: annotation degradation benchmark ─────────────────────────────────────
rule benchmark_degradation:
    """
    注释退化鲁棒性基准：逐步损坏 GFF 注释文件，
    评估 OmniQuant 共识在非模式昆虫场景下的退化行为。
    """
    output:
        tsv = "results/benchmark/annotation_degradation.tsv",
        fig = "results/benchmark/figures/degradation_curve_combined.png"
    params:
        gff = os.path.join(BENCH_REFERENCE_DIR, "Drosophila_melanogaster.BDGP6.54.115.chr.gff3"),
        fc_dea = "results/06.differential_expression/featurecounts/deseq2.60d_vs_1d.csv",
        st_dea = "results/06.differential_expression/stringtie/deseq2.60d_vs_1d.csv",
        sa_dea = "results/06.differential_expression/salmon/deseq2.60d_vs_1d.csv",
        ka_dea = "results/06.differential_expression/kallisto/deseq2.60d_vs_1d.csv",
        fc_counts = "results/05.quantification/matrices/featurecounts/featurecounts_gene_counts_matrix.tsv",
        levels = "0,0.25,0.50,0.75"
    conda:
        "../../envs/dea.yaml"
    benchmark:
        "benchmarks/benchmark_degradation.tsv"
    log:
        "logs/benchmark/degradation.log"
    shell:
        """
        Rscript workflow/scripts/benchmark_annotation_degradation.R \
            --gff {params.gff} \
            --featurecounts-dea {params.fc_dea} \
            --stringtie-dea {params.st_dea} \
            --salmon-dea {params.sa_dea} \
            --kallisto-dea {params.ka_dea} \
            --featurecounts-counts {params.fc_counts} \
            --output-tsv {output.tsv} \
            --output-figure {output.fig} \
            --seeds 5 \
            --levels {params.levels}
        """

# ── Rule: run all benchmarks ───────────────────────────────────────────────────
rule benchmark_all:
    """
    运行所有基准测试 — 子采样稳定性 + 注释退化
    """
    input:
        "results/benchmark/subsampling_stability.tsv",
        "results/benchmark/annotation_degradation.tsv"
