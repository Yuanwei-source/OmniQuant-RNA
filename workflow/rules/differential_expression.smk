# Differential Expression Analysis Rules

DEA_IMPORT_CONFIG = config.get("dea_import", {})
DEA_CONFIG = config.get("dea", {})
DEA_METHODS = [str(method) for method in DEA_CONFIG.get("methods", ["deseq2", "edger", "limma"])]
FEATURECOUNTS_MATRIX = "results/05.quantification/matrices/featurecounts/featurecounts_gene_counts_matrix.tsv"


def build_dea_contrasts():
    configured = DEA_CONFIG.get("comparisons", [])
    if configured == "all":
        sample_table = samples_df if "samples_df" in globals() else pd.read_csv(config["samples"], sep="\t")
        groups = sorted(sample_table["group"].astype(str).drop_duplicates().tolist())
        return [f"{left}_vs_{right}" for idx, left in enumerate(groups) for right in groups[idx + 1:]]
    if isinstance(configured, str):
        return [configured]
    return [str(comp) for comp in configured]


DEA_CONTRASTS = build_dea_contrasts()
DEA_RESULT_TABLES = expand(
    "results/06.differential_expression/{{quantifier}}/{method}.{contrast}.csv",
    method=DEA_METHODS,
    contrast=DEA_CONTRASTS
)

# Helper function to get inputs based on quantifier
def get_dea_primary_input(wildcards):
    if wildcards.quantifier == "featurecounts":
        return FEATURECOUNTS_MATRIX
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
        outdir = directory("results/06.differential_expression/{quantifier}"),
        rds = "results/06.differential_expression/{quantifier}/dea_session.rds",
        norm_counts = "results/06.differential_expression/{quantifier}/normalized_counts.csv",
        import_summary = "results/06.differential_expression/{quantifier}/import_summary.tsv",
        result_tables = DEA_RESULT_TABLES
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
        rds = "results/06.differential_expression/{quantifier}/dea_session.rds"
    output:
        outdir = directory("results/06.differential_expression/{quantifier}/integration"),
        pca = "results/06.differential_expression/{quantifier}/integration/PCA_plot.pdf"
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
        expand("results/06.differential_expression/{quantifier}/integration/PCA_plot.pdf", 
               quantifier=["featurecounts", "stringtie", "salmon", "kallisto"])
