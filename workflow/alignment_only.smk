# Alignment Only Workflow
# Run only read alignment

import pandas as pd

configfile: "config/config.yaml"
samples_df = pd.read_csv(config["samples"], sep="\t")
SAMPLES = samples_df["sample"].tolist()

include: "workflow/rules/alignment.smk"

rule all_alignment:
    input:
        # HISAT2 results
        expand("results/aligned/hisat2/{sample}.bam", sample=SAMPLES),
        # STAR results (optional)
        # expand("results/aligned/star/{sample}Aligned.sortedByCoord.out.bam", sample=SAMPLES)
