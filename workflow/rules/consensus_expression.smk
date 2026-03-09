CONSENSUS_CONFIG = config.get("consensus", {})
CONSENSUS_ENABLED = CONSENSUS_CONFIG.get("enabled", False)
CONSENSUS_QUANTIFIERS = CONSENSUS_CONFIG.get(
    "quantifiers",
    ["featurecounts", "stringtie", "salmon", "kallisto"]
)
CONSENSUS_CONTRASTS = CONSENSUS_CONFIG.get("contrasts", [])
CONSENSUS_METHOD = CONSENSUS_CONFIG.get("methods", ["deseq2"])[0]


def get_consensus_inputs(wildcards):
    return {
        quantifier: f"results/05.differential_expression/{quantifier}/{CONSENSUS_METHOD}.{wildcards.contrast}.csv"
        for quantifier in CONSENSUS_QUANTIFIERS
    }


if CONSENSUS_ENABLED and CONSENSUS_CONTRASTS:
    rule run_consensus_dea:
        """
        Integrate DESeq2 results across quantifiers into a consensus DEA layer.
        """
        input:
            unpack(get_consensus_inputs)
        output:
            table = "results/06.consensus_expression/{contrast}/consensus_results.tsv",
            summary = "results/06.consensus_expression/{contrast}/consensus_summary.tsv",
            scatter = "results/06.consensus_expression/{contrast}/logFC_scatter_salmon_vs_featurecounts.pdf"
        conda:
            "../../envs/dea.yaml"
        log:
            "logs/dea/consensus_{contrast}.log"
        script:
            "../scripts/run_consensus_dea.R"


rule consensus_all:
    input:
        expand(
            "results/06.consensus_expression/{contrast}/consensus_results.tsv",
            contrast=CONSENSUS_CONTRASTS
        )