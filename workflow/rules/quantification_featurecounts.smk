# Quantification Rules
# Gene quantification using featureCounts

rule featurecounts_single:
    """
    Count reads mapped to genes using featureCounts (single sample)
    """
    input:
        bam="results/03.alignment/{sample}.bam",
        gtf=config["reference"]["gtf"]
    output:
        counts="results/04.quantification/featurecounts/{sample}/counts.txt",
        summary="results/04.quantification/featurecounts/{sample}/counts.txt.summary"
    conda:
        "../../envs/featurecounts.yaml"
    params:
        extra=config.get("featurecounts", {}).get("extra", "-p -B -C"),
        feature_type=config.get("featurecounts", {}).get("feature_type", "exon"),
        attribute=config.get("featurecounts", {}).get("attribute", "gene_id")
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

rule featurecounts_all:
    """
    Count reads mapped to genes using featureCounts (all samples together)
    """
    input:
        bams=expand("results/03.alignment/{sample}.bam", sample=SAMPLES),
        gtf=config["reference"]["gtf"]
    output:
        counts="results/04.quantification/featurecounts/all_samples/counts_matrix.txt",
        summary="results/04.quantification/featurecounts/all_samples/counts_matrix.txt.summary"
    conda:
        "../../envs/featurecounts.yaml"
    params:
        extra=config.get("featurecounts", {}).get("extra", "-p -B -C"),
        feature_type=config.get("featurecounts", {}).get("feature_type", "exon"),
        attribute=config.get("featurecounts", {}).get("attribute", "gene_id")
    log:
        "logs/featurecounts/all_samples.log"
    threads: 12
    shell:
        """
        featureCounts \
        -a {input.gtf} \
        -o {output.counts} \
        -t {params.feature_type} \
        -g {params.attribute} \
        -T {threads} \
        {params.extra} \
        {input.bams} 2> {log}
        """

rule featurecounts_transcript:
    """
    Count reads mapped to transcripts using featureCounts
    """
    input:
        bam="results/03.alignment/{sample}.bam",
        gtf=config["reference"]["gtf"]
    output:
        counts="results/04.quantification/featurecounts/{sample}/transcript_counts.txt",
        summary="results/04.quantification/featurecounts/{sample}/transcript_counts.txt.summary"
    conda:
        "../../envs/featurecounts.yaml"
    params:
        extra=config.get("featurecounts", {}).get("extra", "-p -B -C"),
        feature_type="transcript",
        attribute="transcript_id"
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
        bam="results/03.alignment/{sample}.bam",
        gtf=config["reference"]["gtf"]
    output:
        counts="results/04.quantification/featurecounts/{sample}/exon_counts.txt",
        summary="results/04.quantification/featurecounts/{sample}/exon_counts.txt.summary"
    conda:
        "../../envs/featurecounts.yaml"
    params:
        extra=config.get("featurecounts", {}).get("extra", "-p -B -C -f"),
        feature_type="exon",
        attribute="exon_id"
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
        "results/04.quantification/featurecounts/{sample}/counts.txt"
    output:
        "results/quantification_results/{sample}/featurecounts.txt"
    shell:
        "mkdir -p $(dirname {output}) && ln -sf ../../05.quantification/featurecounts/{wildcards.sample}/counts.txt {output}"

rule aggregate_featurecounts_summary:
    """
    Aggregate featureCounts results across all samples
    """
    input:
        expand("results/04.quantification/featurecounts/{sample}/counts.txt", sample=SAMPLES)
    output:
        counts="results/04.quantification/featurecounts/all_samples_counts_matrix.txt"
    conda:
        "../../envs/qc.yaml"
    params:
        samples=SAMPLES,
        input_dir="results/04.quantification/featurecounts"
    log:
        "logs/featurecounts/aggregate_summary.log"
    shell:
        """
        python workflow/scripts/aggregate_featurecounts.py \
            --input-dir {params.input_dir} \
            --samples {params.samples} \
            --output-counts {output.counts} 2> {log}
        """

rule aggregate_featurecounts:
    """
    Create the expected all_samples_counts.txt file (for backward compatibility)
    """
    input:
        "results/04.quantification/featurecounts/all_samples_counts_matrix.txt"
    output:
        "results/04.quantification/featurecounts/all_samples_counts.txt"
    shell:
        "cp {input} {output}"
