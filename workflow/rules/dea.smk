# Differential Expression Analysis Rules

# Helper function to get inputs based on quantifier
def get_dea_input(wildcards):
    if wildcards.quantifier == "featurecounts":
        return "results/04.quantification/featurecounts/all_samples/counts_matrix.txt"
    elif wildcards.quantifier == "stringtie":
        return "results/04.quantification/stringtie/all_samples_gene_counts_matrix.txt"
    elif wildcards.quantifier == "salmon":
        return expand("results/04.quantification/salmon/{sample}/quant.sf", sample=SAMPLES)
    elif wildcards.quantifier == "kallisto":
        return expand("results/04.quantification/kallisto/{sample}/abundance.tsv", sample=SAMPLES)
    else:
        raise ValueError(f"Unknown quantifier: {wildcards.quantifier}")

rule run_dea:
    """
    Run Differential Expression Analysis for a specific quantifier
    """
    input:
        counts = get_dea_input,
        gtf = config["reference"]["gtf"]
    output:
        outdir = directory("results/05.differential_expression/{quantifier}"),
        rds = "results/05.differential_expression/{quantifier}/dea_session.rds",
        norm_counts = "results/05.differential_expression/{quantifier}/normalized_counts.csv"
    conda:
        "../../envs/dea.yaml"
    params:
        sample_file = config["samples"]
    log:
        "logs/dea/{quantifier}_run.log"
    script:
        "../scripts/run_dea.R"

rule integrate_dea:
    """
    Integrate and visualize DEA results for a specific quantifier
    """
    input:
        rds = "results/05.differential_expression/{quantifier}/dea_session.rds"
    output:
        outdir = directory("results/05.differential_expression/{quantifier}/integration"),
        pca = "results/05.differential_expression/{quantifier}/integration/PCA_plot.pdf"
    conda:
        "../../envs/dea.yaml"
    log:
        "logs/dea/{quantifier}_integrate.log"
    script:
        "../scripts/integrate_results.R"

rule dea_all:
    """
    Run DEA for all configured methods
    """
    input:
        expand("results/05.differential_expression/{quantifier}/integration/PCA_plot.pdf", 
               quantifier=["featurecounts", "stringtie", "salmon", "kallisto"])
