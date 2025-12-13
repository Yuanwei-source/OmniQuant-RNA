#!/usr/bin/env python3

"""
Modular workflow for gene quantification using featureCounts
"""

import os
from snakemake.utils import min_version

min_version("7.0")

# Load configuration
configfile: "config/config.yaml"

# Load samples
import pandas as pd
samples = pd.read_csv(config["samples"], sep="\t")["sample"].tolist()
config["samples"] = samples

# Include rules
include: "workflow/rules/quantification_featurecounts.smk"

# Default target rule
rule all:
    input:
        # Individual sample counts
        expand("results/quantification/featurecounts/{sample}/counts.txt", sample=samples),
        expand("results/quantification/featurecounts/{sample}/transcript_counts.txt", sample=samples),
        expand("results/quantification/featurecounts/{sample}/exon_counts.txt", sample=samples),
        # Combined count matrix
        "results/quantification/featurecounts/all_samples/counts_matrix.txt",
        # Compatibility outputs
        expand("results/quantification/{sample}/featurecounts.txt", sample=samples)
