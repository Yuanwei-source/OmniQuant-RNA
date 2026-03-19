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

# Get selected aligner from config
ALIGNER = config.get("aligner", "hisat2")

include: "../rules/decontam.smk"
include: "../rules/quantification_featurecounts.smk"
include: "../rules/quantification_stringtie.smk"
include: "../rules/quantification_kallisto.smk"
include: "../rules/quantification_salmon.smk"

rule all_quantification:
    input:
        # Kallisto results
        expand("results/05.quantification/native/kallisto/per_sample/{sample}/abundance.tsv", sample=SAMPLES),
        # Salmon results  
        expand("results/05.quantification/native/salmon/per_sample/{sample}/quant.sf", sample=SAMPLES),
        # Compatibility symlinks
        expand("results/quantification/{sample}/abundance.tsv", sample=SAMPLES)
