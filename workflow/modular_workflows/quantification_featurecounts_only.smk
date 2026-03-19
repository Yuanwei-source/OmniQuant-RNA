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

# Get selected aligner from config
ALIGNER = config.get("aligner", "hisat2")

include: "../rules/quantification_featurecounts.smk"

rule all:
    input:
        # Individual sample counts
        expand("results/05.quantification/native/featurecounts/per_sample/{sample}/counts.txt", sample=SAMPLES),
        expand("results/05.quantification/native/featurecounts/per_sample/{sample}/transcript_counts.txt", sample=SAMPLES),
        expand("results/05.quantification/native/featurecounts/per_sample/{sample}/exon_counts.txt", sample=SAMPLES),
        # Combined count matrix
        "results/05.quantification/matrices/featurecounts/featurecounts_gene_counts_matrix.tsv",
        # Compatibility outputs
        expand("results/quantification/{sample}/featurecounts.txt", sample=SAMPLES)
