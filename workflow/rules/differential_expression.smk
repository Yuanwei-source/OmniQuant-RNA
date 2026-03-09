# Differential Expression Analysis Rules

DEA_IMPORT_CONFIG = config.get("dea_import", {})

# Helper function to get inputs based on quantifier
def get_dea_primary_input(wildcards):
    if wildcards.quantifier == "featurecounts":
        return "results/04.quantification/featurecounts/all_samples/counts_matrix.txt"
    elif wildcards.quantifier == "stringtie":
        return STRINGTIE_MANIFEST
    elif wildcards.quantifier == "salmon":
        return SALMON_MANIFEST
    elif wildcards.quantifier == "kallisto":
        return KALLISTO_MANIFEST
    else:
        raise ValueError(f"Unknown quantifier: {wildcards.quantifier}")


def get_dea_input_mode(wildcards):
    if wildcards.quantifier == "featurecounts":
        return "gene_counts_matrix"
    return f"tximport_{wildcards.quantifier}"

rule run_dea:
    """
    Run Differential Expression Analysis for a specific quantifier
    """
    wildcard_constraints:
        quantifier = "[^/]+"
    input:
        primary = get_dea_primary_input,
        sample_file = config["samples"],
        tx2gene_master = TX2GENE_MASTER,
        gene_namespace = GENE_NAMESPACE
    output:
        outdir = directory("results/05.differential_expression/{quantifier}"),
        rds = "results/05.differential_expression/{quantifier}/dea_session.rds",
        norm_counts = "results/05.differential_expression/{quantifier}/normalized_counts.csv",
        import_summary = "results/05.differential_expression/{quantifier}/import_summary.tsv"
    conda:
        "../../envs/dea.yaml"
    params:
        input_mode = get_dea_input_mode,
        read_length = config.get("read_length", 150),
        main_only = DEA_IMPORT_CONFIG.get("main_only", True),
        mapping_policy = DEA_IMPORT_CONFIG.get("mapping_policy", "conservative"),
        novel_policy = DEA_IMPORT_CONFIG.get("novel_policy", "discovery_only")
    log:
        "logs/differential_expression/{quantifier}_run.log"
    script:
        "../scripts/perform_quantifier_dea.R"

rule integrate_dea:
    """
    Integrate and visualize DEA results for a specific quantifier
    """
    wildcard_constraints:
        quantifier = "[^/]+"
    input:
        rds = "results/05.differential_expression/{quantifier}/dea_session.rds"
    output:
        outdir = directory("results/05.differential_expression/{quantifier}/integration"),
        pca = "results/05.differential_expression/{quantifier}/integration/PCA_plot.pdf"
    conda:
        "../../envs/dea.yaml"
    log:
        "logs/differential_expression/{quantifier}_integrate.log"
    script:
        "../scripts/integrate_results.R"

rule dea_all:
    """
    Run DEA for all configured methods
    """
    input:
        expand("results/05.differential_expression/{quantifier}/integration/PCA_plot.pdf", 
               quantifier=["featurecounts", "stringtie", "salmon", "kallisto"])
