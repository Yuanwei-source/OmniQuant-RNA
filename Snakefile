# Snakemake workflow for RNA-seq transcriptome quantification
# Author: Yuanwei
# Date: 2025-08-09

include: "workflow/rules/common.smk"

# Include all modular rules
include: "workflow/rules/annotation_conversion.smk"
include: "workflow/rules/qc.smk"
include: "workflow/rules/decontam.smk"
include: "workflow/rules/alignment.smk"
include: "workflow/rules/reference_namespace.smk"
include: "workflow/rules/quantification_kallisto.smk"
include: "workflow/rules/quantification_salmon.smk"
include: "workflow/rules/quantification_featurecounts.smk"
include: "workflow/rules/quantification_stringtie.smk"
include: "workflow/rules/report.smk"
include: "workflow/rules/differential_expression.smk"
include: "workflow/rules/consensus_expression.smk"

# All target outputs - comprehensive RNA-seq pipeline
rule all:
    input:
        # Unified namespace and import manifests
        REFERENCE_TX2GENE,
        GENE_NAMESPACE,
        STRINGTIE_BRIDGE,
        TX2GENE_MASTER,
        SALMON_MANIFEST,
        KALLISTO_MANIFEST,
        STRINGTIE_MANIFEST,

        # Annotation format conversion
        "data/reference/annotation_conversion_complete.flag",
        
        # Transcriptome extraction
        REFERENCE_TRANSCRIPTOME,
        
        # Raw data quality control
        expand("results/01.raw_qc/{sample}_{read}_fastqc.html", sample=SAMPLES, read=["R1", "R2"]),
        expand("results/01.raw_qc/{sample}_{read}_fastqc.zip", sample=SAMPLES, read=["R1", "R2"]),
        
        # Quality trimming with fastp
        expand("results/02.trimmed_data/{sample}_{read}_trimmed.fastq.gz", sample=SAMPLES, read=["R1", "R2"]),
        expand("results/02.trimmed_data/{sample}.fastp.html", sample=SAMPLES),
        
        # Post-trimming quality control
        expand("results/02.trimmed_data/{sample}_{read}_trimmed_fastqc.html", sample=SAMPLES, read=["R1", "R2"]),

        # Optional decontamination outputs
        get_decontam_all_targets(SAMPLES),
        
        # Reference indices
        "data/reference/kallisto_index/transcriptome.idx",
        "data/reference/salmon_index",
        "data/reference/hisat2_index",
        
        # Alignment results
        expand("results/04.alignment/{sample}.bam", sample=SAMPLES),
        expand("results/04.alignment/{sample}.bam.bai", sample=SAMPLES),
        
        # Quantification results - multiple methods
        # Kallisto quantification
        expand("results/05.quantification/native/kallisto/per_sample/{sample}/abundance.tsv", sample=SAMPLES),
        expand("results/05.quantification/native/kallisto/per_sample/{sample}/abundance.h5", sample=SAMPLES),
        "results/05.quantification/matrices/kallisto/kallisto_transcript_counts_matrix.tsv",
        "results/05.quantification/matrices/kallisto/kallisto_transcript_tpm_matrix.tsv",
        
        # Salmon quantification
        expand("results/05.quantification/native/salmon/per_sample/{sample}/quant.sf", sample=SAMPLES),
        "results/05.quantification/matrices/salmon/salmon_transcript_counts_matrix.tsv",
        "results/05.quantification/matrices/salmon/salmon_transcript_tpm_matrix.tsv",
        "results/05.quantification/matrices/salmon/salmon_gene_counts_matrix.tsv",
        "results/05.quantification/matrices/salmon/salmon_gene_tpm_matrix.tsv",
        
        # featureCounts quantification
        "data/reference/genome.featurecounts.gtf",
        "data/reference/genome.featurecounts.summary.tsv",
        expand("results/05.quantification/native/featurecounts/per_sample/{sample}/counts.txt", sample=SAMPLES),
        "results/05.quantification/matrices/featurecounts/featurecounts_gene_counts_matrix.tsv",
        
        # StringTie quantification
        expand("results/05.quantification/native/stringtie/per_sample/{sample}/final/transcripts.gtf", sample=SAMPLES),
        expand("results/05.quantification/native/stringtie/per_sample/{sample}/final/gene_abundances.tab", sample=SAMPLES),
        "results/05.quantification/native/stringtie/merged/merged.gtf",
        "results/05.quantification/audit/stringtie/gene_id_mapping.tsv",
        "results/05.quantification/matrices/stringtie/stringtie_gene_counts_matrix.tsv",
        "results/05.quantification/matrices/stringtie/stringtie_transcript_counts_matrix.tsv",
        "results/05.quantification/matrices/stringtie/stringtie_gene_tpm_matrix.tsv",
        "results/05.quantification/matrices/stringtie/stringtie_transcript_tpm_matrix.tsv",
        "results/05.quantification/matrices/stringtie/stringtie_gene_fpkm_matrix.tsv",
        "results/05.quantification/matrices/stringtie/stringtie_transcript_fpkm_matrix.tsv",
        
        # StringTie results with original gene IDs (for downstream analysis)
        "results/05.quantification/matrices/stringtie/stringtie_original_gene_counts_matrix.tsv",
        "results/05.quantification/matrices/stringtie/stringtie_original_transcript_counts_matrix.tsv",
        "results/05.quantification/matrices/stringtie/stringtie_original_gene_tpm_matrix.tsv",
        "results/05.quantification/matrices/stringtie/stringtie_original_transcript_tpm_matrix.tsv",
        
        # DEA results
        expand("results/06.differential_expression/{quantifier}/integration/PCA_plot.pdf", quantifier=["featurecounts", "stringtie", "salmon", "kallisto"]),

        # Consensus DEA results
        expand("results/07.consensus_expression/{contrast}/consensus_results.tsv", contrast=CONSENSUS_CONTRASTS),
        expand("results/07.consensus_expression/{contrast}/consensus_summary.tsv", contrast=CONSENSUS_CONTRASTS),
        expand("results/07.consensus_expression/{contrast}/tier_diagnostics.tsv", contrast=CONSENSUS_CONTRASTS),
        expand("results/07.consensus_expression/{contrast}/significance_membership.tsv", contrast=CONSENSUS_CONTRASTS),
        expand("results/07.consensus_expression/{contrast}/sensitivity_analysis.tsv", contrast=CONSENSUS_CONTRASTS),
        expand("results/07.consensus_expression/{contrast}/logFC_scatter_salmon_vs_featurecounts.pdf", contrast=CONSENSUS_CONTRASTS),
        expand("results/07.consensus_expression/{contrast}/consensus_volcano.pdf", contrast=CONSENSUS_CONTRASTS),
        expand("results/07.consensus_expression/{contrast}/significance_upset.pdf", contrast=CONSENSUS_CONTRASTS),
        
        # Reports and summaries
        "results/08.reports/multiqc_report.html"

# Optional rules for different analysis subsets

rule qc_only:
    """Run only quality control steps"""
    input:
        expand("results/01.raw_qc/{sample}_{read}_fastqc.html", sample=SAMPLES, read=["R1", "R2"]),
        expand("results/02.trimmed_data/{sample}_{read}_trimmed.fastq.gz", sample=SAMPLES, read=["R1", "R2"]),
        expand("results/02.trimmed_data/{sample}_{read}_trimmed_fastqc.html", sample=SAMPLES, read=["R1", "R2"]),
        get_decontam_all_targets(SAMPLES)

rule stringtie_only:
    """Run only StringTie quantification with original gene IDs"""
    input:
        expand("results/05.quantification/native/stringtie/per_sample/{sample}/final/transcripts.gtf", sample=SAMPLES),
        expand("results/05.quantification/native/stringtie/per_sample/{sample}/final/gene_abundances.tab", sample=SAMPLES),
        "results/05.quantification/native/stringtie/merged/merged.gtf",
        "results/05.quantification/audit/stringtie/gene_id_mapping.tsv",
        "results/05.quantification/matrices/stringtie/stringtie_original_gene_counts_matrix.tsv",
        "results/05.quantification/matrices/stringtie/stringtie_original_gene_tpm_matrix.tsv"

rule alignment_only:
    """Run quality control and alignment only"""
    input:
        get_decontam_all_targets(SAMPLES),
        expand("results/04.alignment/{sample}.bam", sample=SAMPLES),
        expand("results/04.alignment/{sample}.bam.bai", sample=SAMPLES)

rule quantification_kallisto_only:
    """Run Kallisto quantification only"""
    input:
        "data/reference/kallisto_index/transcriptome.idx",
        expand("results/05.quantification/native/kallisto/per_sample/{sample}/abundance.tsv", sample=SAMPLES),
        "results/05.quantification/matrices/kallisto/kallisto_transcript_counts_matrix.tsv",
        "results/05.quantification/matrices/kallisto/kallisto_transcript_tpm_matrix.tsv"

rule quantification_salmon_only:
    """Run Salmon quantification only"""
    input:
        "data/reference/salmon_index",
        expand("results/05.quantification/native/salmon/per_sample/{sample}/quant.sf", sample=SAMPLES),
        "results/05.quantification/matrices/salmon/salmon_transcript_counts_matrix.tsv",
        "results/05.quantification/matrices/salmon/salmon_transcript_tpm_matrix.tsv",
        "results/05.quantification/matrices/salmon/salmon_gene_counts_matrix.tsv",
        "results/05.quantification/matrices/salmon/salmon_gene_tpm_matrix.tsv"

rule consensus_only:
    """Run configured consensus DEA contrasts only"""
    input:
        expand("results/07.consensus_expression/{contrast}/consensus_results.tsv", contrast=CONSENSUS_CONTRASTS)

rule quantification_featurecounts_only:
    """Run featureCounts quantification only"""
    input:
        "data/reference/genome.featurecounts.gtf",
        "data/reference/genome.featurecounts.summary.tsv",
        expand("results/05.quantification/native/featurecounts/per_sample/{sample}/counts.txt", sample=SAMPLES),
        "results/05.quantification/matrices/featurecounts/featurecounts_gene_counts_matrix.tsv"

rule quantification_stringtie_only:
    """Run StringTie quantification only"""
    input:
        expand("results/05.quantification/native/stringtie/per_sample/{sample}/final/transcripts.gtf", sample=SAMPLES),
        expand("results/05.quantification/native/stringtie/per_sample/{sample}/final/gene_abundances.tab", sample=SAMPLES),
        "results/05.quantification/native/stringtie/merged/merged.gtf",
        "results/05.quantification/matrices/stringtie/stringtie_gene_counts_matrix.tsv",
        "results/05.quantification/matrices/stringtie/stringtie_transcript_counts_matrix.tsv",
        "results/05.quantification/matrices/stringtie/stringtie_gene_tpm_matrix.tsv",
        "results/05.quantification/matrices/stringtie/stringtie_transcript_tpm_matrix.tsv",
        "results/05.quantification/matrices/stringtie/stringtie_gene_fpkm_matrix.tsv",
        "results/05.quantification/matrices/stringtie/stringtie_transcript_fpkm_matrix.tsv"
