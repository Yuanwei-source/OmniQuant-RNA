import os

FEATURECOUNTS_CONFIG = config.get("featurecounts", {})
FEATURECOUNTS_GTF_SUMMARY = f"{REFERENCE_DIR}/genome.featurecounts.summary.tsv"

def resolve_featurecounts_attribute(config_value, target_type="gene"):
    if config_value and str(config_value).strip().lower() != "auto":
        return config_value

    defaults = {
        "gene": "gene_id",
        "transcript": "transcript_id",
        "exon": "exon_id"
    }
    return defaults.get(target_type, "gene_id")

# Quantification Rules
# Gene quantification using featureCounts

rule normalize_featurecounts_annotation:
    """
    Normalize GTF attributes for featureCounts compatibility
    """
    input:
        gtf=REFERENCE_GTF
    output:
        gtf=FEATURECOUNTS_GTF,
        summary=FEATURECOUNTS_GTF_SUMMARY
    conda:
        "../../envs/qc.yaml"
    threads: 8
    log:
        "logs/featurecounts/normalize_annotation.log"
    shell:
        """
        python3 workflow/scripts/normalize_featurecounts_gtf.py \
            --input {input.gtf} \
            --output {output.gtf} \
            --summary {output.summary} \
            --threads {threads} > {log} 2>&1
        """

rule featurecounts_single:
    """
    Count reads mapped to genes using featureCounts (single sample)
    """
    input:
        bam="results/04.alignment/{sample}.bam",
        gtf=FEATURECOUNTS_GTF,
        summary=FEATURECOUNTS_GTF_SUMMARY
    output:
        counts="results/05.quantification/native/featurecounts/per_sample/{sample}/counts.txt",
        summary="results/05.quantification/native/featurecounts/per_sample/{sample}/counts.txt.summary"
    conda:
        "../../envs/featurecounts.yaml"
    params:
        extra=FEATURECOUNTS_CONFIG.get("extra", "-p -B -C"),
        feature_type=FEATURECOUNTS_CONFIG.get("feature_type", "exon"),
        attribute=resolve_featurecounts_attribute(FEATURECOUNTS_CONFIG.get("attribute", "auto"), "gene")
    log:
        "logs/featurecounts/{sample}.log"
    threads: 8
    shell:
        """
        featureCounts \
        -a {input.gtf} \
        -o {output.counts} \
        -t {params.feature_type} \
        -g {params.attribute} \
        -T {threads} \
        {params.extra} \
        {input.bam} 2> {log}
        """

rule featurecounts_transcript:
    """
    Count reads mapped to transcripts using featureCounts
    """
    input:
        bam="results/04.alignment/{sample}.bam",
        gtf=FEATURECOUNTS_GTF,
        summary=FEATURECOUNTS_GTF_SUMMARY
    output:
        counts="results/05.quantification/native/featurecounts/per_sample/{sample}/transcript_counts.txt",
        summary="results/05.quantification/native/featurecounts/per_sample/{sample}/transcript_counts.txt.summary"
    conda:
        "../../envs/featurecounts.yaml"
    params:
        extra=FEATURECOUNTS_CONFIG.get("extra", "-p -B -C"),
        feature_type="transcript",
        attribute=resolve_featurecounts_attribute(FEATURECOUNTS_CONFIG.get("transcript_attribute", "auto"), "transcript")
    log:
        "logs/featurecounts/{sample}_transcript.log"
    threads: 8
    shell:
        """
        featureCounts \
        -a {input.gtf} \
        -o {output.counts} \
        -t {params.feature_type} \
        -g {params.attribute} \
        -T {threads} \
        {params.extra} \
        {input.bam} 2> {log}
        """

rule featurecounts_exon:
    """
    Count reads mapped to exons using featureCounts
    """
    input:
        bam="results/04.alignment/{sample}.bam",
        gtf=FEATURECOUNTS_GTF,
        summary=FEATURECOUNTS_GTF_SUMMARY
    output:
        counts="results/05.quantification/native/featurecounts/per_sample/{sample}/exon_counts.txt",
        summary="results/05.quantification/native/featurecounts/per_sample/{sample}/exon_counts.txt.summary"
    conda:
        "../../envs/featurecounts.yaml"
    params:
        extra=FEATURECOUNTS_CONFIG.get("extra", "-p -B -C -f"),
        feature_type="exon",
        attribute=resolve_featurecounts_attribute(FEATURECOUNTS_CONFIG.get("exon_attribute", "auto"), "exon")
    log:
        "logs/featurecounts/{sample}_exon.log"
    threads: 8
    shell:
        """
        featureCounts \
        -a {input.gtf} \
        -o {output.counts} \
        -t {params.feature_type} \
        -g {params.attribute} \
        -T {threads} \
        {params.extra} \
        {input.bam} 2> {log}
        """

# Compatibility rule for main workflow
rule quantification_results_featurecounts:
    """
    Create symlink for main workflow compatibility
    """
    input:
        "results/05.quantification/native/featurecounts/per_sample/{sample}/counts.txt"
    output:
        "results/quantification_results/{sample}/featurecounts.txt"
    shell:
        "mkdir -p $(dirname {output}) && ln -sf ../../05.quantification/native/featurecounts/per_sample/{wildcards.sample}/counts.txt {output}"

rule aggregate_featurecounts_summary:
    """
    Aggregate featureCounts per-sample results into the canonical gene count matrix
    """
    input:
        expand("results/05.quantification/native/featurecounts/per_sample/{sample}/counts.txt", sample=SAMPLES)
    output:
        counts="results/05.quantification/matrices/featurecounts/featurecounts_gene_counts_matrix.tsv"
    conda:
        "../../envs/qc.yaml"
    params:
        samples=SAMPLES,
        input_dir="results/05.quantification/native/featurecounts/per_sample"
    log:
        "logs/featurecounts/aggregate_summary.log"
    shell:
        """
        python3 workflow/scripts/aggregate_featurecounts.py \
            --input-dir {params.input_dir} \
            --samples {params.samples} \
            --output-counts {output.counts} 2> {log}
        """
