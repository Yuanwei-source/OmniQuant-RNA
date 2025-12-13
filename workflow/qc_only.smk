# Quality Control Only Workflow
# Run only FastQC analysis

import pandas as pd

configfile: "config/config.yaml"
samples_df = pd.read_csv(config["samples"], sep="\t")
SAMPLES = samples_df["sample"].tolist()

include: "workflow/rules/qc.smk"

rule all_qc:
    input:
        expand("results/fastqc/{sample}_{read}_fastqc.html", sample=SAMPLES, read=["R1", "R2"]),
        expand("results/fastqc/{sample}_{read}_fastqc.zip", sample=SAMPLES, read=["R1", "R2"])
