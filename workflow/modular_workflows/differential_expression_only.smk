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

include: "../rules/reference_namespace.smk"
include: "../rules/differential_expression.smk"

rule all_differential_expression:
    input:
        rules.dea_all.input
