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

# Functions to get input files from samples.tsv
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

# Get selected aligner from config
ALIGNER = config.get("aligner", "hisat2")

include: "../rules/qc.smk"

rule all_qc:
    input:
        expand("results/fastqc/{sample}_{read}_fastqc.html", sample=SAMPLES, read=["R1", "R2"]),
        expand("results/fastqc/{sample}_{read}_fastqc.zip", sample=SAMPLES, read=["R1", "R2"])
