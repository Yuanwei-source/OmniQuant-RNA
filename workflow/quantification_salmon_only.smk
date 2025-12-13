# Salmon Quantification Only Workflow
# Run only Salmon transcript quantification

import pandas as pd

configfile: "config/config.yaml"
samples_df = pd.read_csv(config["samples"], sep="\t")
SAMPLES = samples_df["sample"].tolist()

include: "rules/quantification_salmon.smk"

rule all_salmon:
    input:
        # Salmon results  
        expand("results/quantification/salmon/{sample}/quant.sf", sample=SAMPLES),
        # Compatibility symlinks
        expand("results/quantification/{sample}/quant.sf", sample=SAMPLES)
