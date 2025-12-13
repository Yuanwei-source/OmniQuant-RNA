#!/usr/bin/env python3

"""
Modular workflow for transcript quantification using StringTie
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
include: "workflow/rules/quantification_stringtie.smk"

# Default target rule
rule all:
    input:
        # StringTie outputs
        expand("results/quantification/stringtie/{sample}/transcripts.gtf", sample=samples),
        expand("results/quantification/stringtie/{sample}/gene_abundances.tab", sample=samples),
        "results/quantification/stringtie/merged.gtf",
        expand("results/quantification/stringtie/{sample}/final_transcripts.gtf", sample=samples),
        expand("results/quantification/stringtie/{sample}/final_gene_abundances.tab", sample=samples),
        # Compatibility outputs
        expand("results/quantification/{sample}/stringtie_abundances.tab", sample=samples)
