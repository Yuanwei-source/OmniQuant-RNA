# Kallisto Quantification Only Workflow
# Run only Kallisto transcript quantification

import pandas as pd

configfile: "config/config.yaml"
samples_df = pd.read_csv(config["samples"], sep="\t")
SAMPLES = samples_df["sample"].tolist()

include: "rules/quantification_kallisto.smk"

rule all_kallisto:
    input:
        # Kallisto results
        expand("results/quantification/kallisto/{sample}/abundance.tsv", sample=SAMPLES),
        # Compatibility symlinks
        expand("results/quantification/{sample}/abundance.tsv", sample=SAMPLES)
