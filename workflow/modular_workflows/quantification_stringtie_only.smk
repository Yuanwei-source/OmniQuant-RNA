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

include: "../rules/quantification_stringtie.smk"

rule all:
    input:
        # StringTie outputs
        expand("results/05.quantification/native/stringtie/per_sample/{sample}/assembly/transcripts.gtf", sample=SAMPLES),
        expand("results/05.quantification/native/stringtie/per_sample/{sample}/assembly/gene_abundances.tab", sample=SAMPLES),
        "results/05.quantification/native/stringtie/merged/merged.gtf",
        expand("results/05.quantification/native/stringtie/per_sample/{sample}/final/transcripts.gtf", sample=SAMPLES),
        expand("results/05.quantification/native/stringtie/per_sample/{sample}/final/gene_abundances.tab", sample=SAMPLES),
        # Compatibility outputs
        expand("results/quantification/{sample}/stringtie_abundances.tab", sample=SAMPLES)
