# Quantification Only Workflow
# Run only transcript quantification

import pandas as pd

configfile: "config/config.yaml"
samples_df = pd.read_csv(config["samples"], sep="\t")
SAMPLES = samples_df["sample"].tolist()

include: "rules/quantification_kallisto.smk"
include: "rules/quantification_salmon.smk"

rule all_quantification:
    input:
        # Kallisto results
        expand("results/quantification/kallisto/{sample}/abundance.tsv", sample=SAMPLES),
        # Salmon results  
        expand("results/quantification/salmon/{sample}/quant.sf", sample=SAMPLES),
        # Compatibility symlinks
        expand("results/quantification/{sample}/abundance.tsv", sample=SAMPLES)
