CONSENSUS_CONFIG = config.get("consensus", {})
CONSENSUS_ENABLED = CONSENSUS_CONFIG.get("enabled", False)
CONSENSUS_QUANTIFIERS = CONSENSUS_CONFIG.get(
    "quantifiers",
    ["featurecounts", "stringtie", "salmon", "kallisto"]
)
CONSENSUS_METHOD = CONSENSUS_CONFIG.get("methods", ["deseq2"])[0]


def build_consensus_contrasts():
    configured = CONSENSUS_CONFIG.get("contrasts", [])
    if configured == "all":
        dea_config = config.get("dea", {})
        dea_comparisons = dea_config.get("comparisons", [])

        if isinstance(dea_comparisons, str) and dea_comparisons != "all":
            return [dea_comparisons]

        if isinstance(dea_comparisons, (list, tuple)):
            explicit = [str(comp) for comp in dea_comparisons if str(comp).strip()]
            if explicit and explicit != ["all"]:
                return explicit

        sample_table = samples_df if "samples_df" in globals() else pd.read_csv(config["samples"], sep="\t")
        groups = sorted(sample_table["group"].astype(str).drop_duplicates().tolist())
        return [f"{left}_vs_{right}" for idx, left in enumerate(groups) for right in groups[idx + 1:]]
    if isinstance(configured, str):
        return [configured]
    return list(configured)


CONSENSUS_CONTRASTS = build_consensus_contrasts()


def get_consensus_inputs(wildcards):
    return {
        quantifier: f"results/06.differential_expression/{quantifier}/{CONSENSUS_METHOD}.{wildcards.contrast}.csv"
        for quantifier in CONSENSUS_QUANTIFIERS
    }


if CONSENSUS_ENABLED and CONSENSUS_CONTRASTS:
    rule run_consensus_dea:
        """
        Integrate DESeq2 results across quantifiers into a consensus DEA layer.
        """
        input:
            unpack(get_consensus_inputs),
            tx2gene_master = "results/00.reference/tx2gene_master.tsv"
        output:
            table = "results/07.consensus_expression/{contrast}/consensus_results.tsv",
            summary = "results/07.consensus_expression/{contrast}/consensus_summary.tsv",
            diagnostics = "results/07.consensus_expression/{contrast}/tier_diagnostics.tsv",
            membership = "results/07.consensus_expression/{contrast}/significance_membership.tsv",
            sensitivity = "results/07.consensus_expression/{contrast}/sensitivity_analysis.tsv",
            scatter = "results/07.consensus_expression/{contrast}/logFC_scatter_salmon_vs_featurecounts.pdf",
            volcano = "results/07.consensus_expression/{contrast}/consensus_volcano.pdf",
            upset = "results/07.consensus_expression/{contrast}/significance_upset.pdf"
        conda:
            "../../envs/dea.yaml"
        log:
            "logs/dea/consensus_{contrast}.log"
        script:
            "../scripts/run_consensus_dea.R"


rule consensus_all:
    input:
        expand(
            "results/07.consensus_expression/{contrast}/consensus_results.tsv",
            contrast=CONSENSUS_CONTRASTS
        )