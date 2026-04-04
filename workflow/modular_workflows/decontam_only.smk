import pandas as pd
import os
import subprocess
import sys

# Configuration file
configfile: "config/config.yaml"

include: "../rules/reference_config.smk"

# Sample information
samples_df = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
SAMPLES = samples_df["sample"].tolist()


def get_r1(wildcards):
    return samples_df.loc[wildcards.sample, "fq1"]


def get_r2(wildcards):
    return samples_df.loc[wildcards.sample, "fq2"]


FASTQ_SUFFIXES = [
    ".fastq.gz",
    ".fq.gz",
    ".fastq.bz2",
    ".fq.bz2",
    ".fastq",
    ".fq",
]


def get_fastq_suffix(path):
    lower_path = str(path).lower()
    for suffix in FASTQ_SUFFIXES:
        if lower_path.endswith(suffix):
            return suffix
    raise ValueError(
        f"Unsupported FASTQ extension for '{path}'. Supported suffixes: {', '.join(FASTQ_SUFFIXES)}"
    )


def get_fastqc_alias_path(sample, read_label, source_path):
    suffix = get_fastq_suffix(source_path)
    return os.path.join(".snakemake", "fastqc_aliases", sample, f"{sample}_{read_label}{suffix}")


def decontam_enabled():
    return bool(config.get("decontam", {}).get("enabled", False))


def get_analysis_reads(wildcards):
    if decontam_enabled():
        return {
            "r1": f"results/03.decontam/clean/{wildcards.sample}_R1_clean.fastq.gz",
            "r2": f"results/03.decontam/clean/{wildcards.sample}_R2_clean.fastq.gz",
        }

    return {
        "r1": f"results/02.trimmed_data/{wildcards.sample}_R1_trimmed.fastq.gz",
        "r2": f"results/02.trimmed_data/{wildcards.sample}_R2_trimmed.fastq.gz",
    }


def get_analysis_r1(wildcards):
    return get_analysis_reads(wildcards)["r1"]


def get_analysis_r2(wildcards):
    return get_analysis_reads(wildcards)["r2"]


include: "../rules/qc.smk"
include: "../rules/decontam.smk"


rule all_decontam:
    input:
        expand("results/03.decontam/clean/{sample}_R1_clean.fastq.gz", sample=SAMPLES),
        expand("results/03.decontam/clean/{sample}_R2_clean.fastq.gz", sample=SAMPLES),
        expand("results/03.decontam/stats/{sample}_decision_summary.tsv", sample=SAMPLES),
        expand("results/03.decontam/qc/{sample}_{read}_clean_fastqc.html", sample=SAMPLES, read=["R1", "R2"]),
        "results/03.decontam/stats/project_decontam_summary.tsv",
        "results/03.decontam/clues/tables/sample_microbial_burden.tsv",
        "results/03.decontam/clues/tables/priority_targets.tsv",
        "results/03.decontam/clues/plots/microbial_composition_stacked_bar.pdf"