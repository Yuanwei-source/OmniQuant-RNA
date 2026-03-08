import os
import gzip

def detect_gtf_attribute(gtf_path, target_type="gene"):
    if not os.path.exists(gtf_path):
        return "gene_id" if target_type == "gene" else "transcript_id"
    
    # Check what attributes exist in the GTF
    attrs = {"gene_id": 0, "gene_name": 0, "gene": 0, "Parent": 0, "transcript_id": 0, "ID": 0}
    open_func = gzip.open if gtf_path.endswith(".gz") else open
    
    try:
        with open_func(gtf_path, "rt") as f:
            lines_checked = 0
            for line in f:
                if line.startswith("#"): continue
                parts = line.strip().split("\t")
                if len(parts) < 9: continue
                
                # For featureCounts -g, we usually look at exon rows for gene_id
                if parts[2] != "exon" and parts[2] != "transcript" and parts[2] != "gene": continue
                
                attr_str = parts[8]
                for k in attrs.keys():
                    if f"{k} " in attr_str or f"{k}=" in attr_str:
                        attrs[k] += 1
                
                lines_checked += 1
                if lines_checked >= 500:
                    break
    except Exception:
        pass
    
    if target_type == "gene":
        if attrs["gene_id"] > 0: return "gene_id"
        if attrs["gene_name"] > 0: return "gene_name"
        if attrs["Parent"] > 0: return "Parent"
        if attrs["gene"] > 0: return "gene"
        return "gene_id"
    else:
        if attrs["transcript_id"] > 0: return "transcript_id"
        if attrs["ID"] > 0: return "ID"
        if attrs["Parent"] > 0: return "Parent"
        return "transcript_id"

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
        attribute=detect_gtf_attribute(config["reference"]["gtf"], "gene")
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
        attribute=detect_gtf_attribute(config["reference"]["gtf"], "gene")
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
        attribute=detect_gtf_attribute(config["reference"]["gtf"], "transcript")
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
