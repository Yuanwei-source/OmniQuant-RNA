import os
import gzip

FEATURECOUNTS_CONFIG = config.get("featurecounts", {})
FEATURECOUNTS_GTF = FEATURECOUNTS_CONFIG.get("normalized_gtf", "data/reference/genome.featurecounts.gtf")
FEATURECOUNTS_GTF_SUMMARY = FEATURECOUNTS_CONFIG.get(
    "normalized_gtf_summary",
    "data/reference/genome.featurecounts.summary.tsv"
)

def parse_annotation_attributes(attr_str):
    attrs = {}
    for raw_field in attr_str.strip().strip(";").split(";"):
        field = raw_field.strip()
        if not field:
            continue

        if "=" in field:
            key, value = field.split("=", 1)
        elif " " in field:
            key, value = field.split(" ", 1)
        else:
            continue

        attrs[key.strip()] = value.strip().strip('"')

    return attrs

def detect_gtf_attribute(gtf_path, target_type="gene"):
    candidate_map = {
        "gene": ("gene_id", "gene_name", "gene", "Parent", "ID"),
        "transcript": ("transcript_id", "ID", "Parent"),
        "exon": ("exon_id", "ID", "Parent", "gene_id", "transcript_id"),
    }
    feature_map = {
        "gene": {"exon", "transcript", "gene"},
        "transcript": {"transcript", "exon"},
        "exon": {"exon"},
    }
    candidates = candidate_map.get(target_type, candidate_map["gene"])
    valid_features = feature_map.get(target_type, feature_map["gene"])

    if not os.path.exists(gtf_path):
        return candidates[0]

    attrs = {key: 0 for key in {attr for values in candidate_map.values() for attr in values}}
    open_func = gzip.open if gtf_path.endswith(".gz") else open

    try:
        with open_func(gtf_path, "rt") as f:
            lines_checked = 0
            for line in f:
                if line.startswith("#"):
                    continue

                parts = line.strip().split("\t")
                if len(parts) < 9:
                    continue
                if parts[2] not in valid_features:
                    continue

                parsed_attrs = parse_annotation_attributes(parts[8])
                for key in candidates:
                    if key in parsed_attrs and parsed_attrs[key]:
                        attrs[key] += 1

                lines_checked += 1
                if lines_checked >= 500:
                    break
    except Exception:
        pass

    for key in candidates:
        if attrs.get(key, 0) > 0:
            return key

    return candidates[0]

def resolve_featurecounts_attribute(config_value, gtf_path, target_type):
    if config_value and str(config_value).strip().lower() != "auto":
        return config_value

    return detect_gtf_attribute(gtf_path, target_type)

GTF_PATH = config.get("reference", {}).get("gtf", "")

FC_GTF_PATH = FEATURECOUNTS_GTF if os.path.exists(FEATURECOUNTS_GTF) else GTF_PATH

FC_GENE = resolve_featurecounts_attribute(
    FEATURECOUNTS_CONFIG.get("attribute", "auto"),
    FC_GTF_PATH,
    "gene"
)
FC_TRANSCRIPT = resolve_featurecounts_attribute(
    FEATURECOUNTS_CONFIG.get("transcript_attribute", "auto"),
    FC_GTF_PATH,
    "transcript"
)
FC_EXON = resolve_featurecounts_attribute(
    FEATURECOUNTS_CONFIG.get("exon_attribute", "auto"),
    FC_GTF_PATH,
    "exon"
)

# Quantification Rules
# Gene quantification using featureCounts

rule normalize_featurecounts_annotation:
    """
    Normalize GTF attributes for featureCounts compatibility
    """
    input:
        gtf=config["reference"]["gtf"]
    output:
        gtf=FEATURECOUNTS_GTF,
        summary=FEATURECOUNTS_GTF_SUMMARY
    conda:
        "../../envs/qc.yaml"
    log:
        "logs/featurecounts/normalize_annotation.log"
    shell:
        """
        python workflow/scripts/normalize_featurecounts_gtf.py \
            --input {input.gtf} \
            --output {output.gtf} \
            --summary {output.summary} > {log} 2>&1
        """

rule featurecounts_single:
    """
    Count reads mapped to genes using featureCounts (single sample)
    """
    input:
        bam="results/03.alignment/{sample}.bam",
        gtf=FEATURECOUNTS_GTF,
        summary=FEATURECOUNTS_GTF_SUMMARY
    output:
        counts="results/04.quantification/featurecounts/{sample}/counts.txt",
        summary="results/04.quantification/featurecounts/{sample}/counts.txt.summary"
    conda:
        "../../envs/featurecounts.yaml"
    params:
        extra=FEATURECOUNTS_CONFIG.get("extra", "-p -B -C"),
        feature_type=FEATURECOUNTS_CONFIG.get("feature_type", "exon"),
        attribute=FC_GENE
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
        gtf=FEATURECOUNTS_GTF,
        summary=FEATURECOUNTS_GTF_SUMMARY
    output:
        counts="results/04.quantification/featurecounts/all_samples/counts_matrix.txt",
        summary="results/04.quantification/featurecounts/all_samples/counts_matrix.txt.summary"
    conda:
        "../../envs/featurecounts.yaml"
    params:
        extra=FEATURECOUNTS_CONFIG.get("extra", "-p -B -C"),
        feature_type=FEATURECOUNTS_CONFIG.get("feature_type", "exon"),
        attribute=FC_GENE
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
        gtf=FEATURECOUNTS_GTF,
        summary=FEATURECOUNTS_GTF_SUMMARY
    output:
        counts="results/04.quantification/featurecounts/{sample}/transcript_counts.txt",
        summary="results/04.quantification/featurecounts/{sample}/transcript_counts.txt.summary"
    conda:
        "../../envs/featurecounts.yaml"
    params:
        extra=FEATURECOUNTS_CONFIG.get("extra", "-p -B -C"),
        feature_type="transcript",
        attribute=FC_TRANSCRIPT
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
        gtf=FEATURECOUNTS_GTF,
        summary=FEATURECOUNTS_GTF_SUMMARY
    output:
        counts="results/04.quantification/featurecounts/{sample}/exon_counts.txt",
        summary="results/04.quantification/featurecounts/{sample}/exon_counts.txt.summary"
    conda:
        "../../envs/featurecounts.yaml"
    params:
        extra=FEATURECOUNTS_CONFIG.get("extra", "-p -B -C -f"),
        feature_type="exon",
        attribute=FC_EXON
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
