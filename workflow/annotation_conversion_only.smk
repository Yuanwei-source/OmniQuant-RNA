# Snakemake workflow for annotation format conversion only
# This workflow converts between GFF3 and GTF annotation formats

import pandas as pd
from snakemake.utils import validate

# Configuration file
configfile: "config/config.yaml"

# Include annotation conversion rules
include: "workflow/rules/annotation_conversion.smk"

# Target output for conversion only
rule all:
    input:
        # Annotation format conversion completion flag
        "data/reference/annotation_conversion_complete.flag"

# Rule to run only annotation conversion
rule conversion_only:
    input:
        "data/reference/annotation_conversion_complete.flag"
    output:
        "results/annotation_conversion_summary.txt"
    shell:
        """
        echo "Annotation format conversion completed successfully" > {output}
        echo "Conversion flag file: {input}" >> {output}
        echo "Date: $(date)" >> {output}
        """
