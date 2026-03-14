# Unified reference namespace and tx2gene rules

NAMESPACE_CONFIG = config.get("namespace", {})
REFERENCE_OUTPUT_DIR = NAMESPACE_CONFIG.get("output_dir", "results/00.reference")
REFERENCE_TX2GENE = NAMESPACE_CONFIG.get("tx2gene_reference", f"{REFERENCE_OUTPUT_DIR}/tx2gene_reference.tsv")
GENE_NAMESPACE = NAMESPACE_CONFIG.get("gene_namespace", f"{REFERENCE_OUTPUT_DIR}/gene_namespace.tsv")
STRINGTIE_BRIDGE = NAMESPACE_CONFIG.get("stringtie_bridge", f"{REFERENCE_OUTPUT_DIR}/stringtie_tx2gene_bridge.tsv")
TX2GENE_MASTER = NAMESPACE_CONFIG.get("tx2gene_master", f"{REFERENCE_OUTPUT_DIR}/tx2gene_master.tsv")
IMPORT_MANIFEST_DIR = NAMESPACE_CONFIG.get("import_manifest_dir", f"{REFERENCE_OUTPUT_DIR}/import_manifests")
SALMON_MANIFEST = f"{IMPORT_MANIFEST_DIR}/salmon.tsv"
KALLISTO_MANIFEST = f"{IMPORT_MANIFEST_DIR}/kallisto.tsv"
STRINGTIE_MANIFEST = f"{IMPORT_MANIFEST_DIR}/stringtie.tsv"


rule build_reference_tx2gene:
    """Build reference transcript-to-gene mapping from the reference GTF."""
    input:
        gtf=config["reference"]["gtf"]
    output:
        REFERENCE_TX2GENE
    conda:
        "../../envs/qc.yaml"
    log:
        "logs/reference_namespace/build_reference_tx2gene.log"
    shell:
        """
        python workflow/scripts/build_reference_tx2gene.py \
            --gtf {input.gtf} \
            --output {output} > {log} 2>&1
        """


rule build_gene_namespace:
    """Build unified gene namespace from the reference tx2gene mapping."""
    input:
        REFERENCE_TX2GENE
    output:
        GENE_NAMESPACE
    conda:
        "../../envs/qc.yaml"
    log:
        "logs/reference_namespace/build_gene_namespace.log"
    shell:
        """
        python workflow/scripts/build_gene_namespace.py \
            --tx2gene {input} \
            --output {output} > {log} 2>&1
        """


rule build_stringtie_bridge:
    """Build conservative StringTie-to-reference mapping bridge."""
    input:
        merged_gtf="results/05.quantification/native/stringtie/merged/merged.gtf",
        reference_tx2gene=REFERENCE_TX2GENE
    output:
        STRINGTIE_BRIDGE
    conda:
        "../../envs/qc.yaml"
    log:
        "logs/reference_namespace/build_stringtie_bridge.log"
    shell:
        """
        python workflow/scripts/build_stringtie_bridge.py \
            --merged-gtf {input.merged_gtf} \
            --reference-tx2gene {input.reference_tx2gene} \
            --output {output} > {log} 2>&1
        """


rule build_tx2gene_master:
    """Combine reference and StringTie mappings into a master tx2gene table."""
    input:
        reference_tx2gene=REFERENCE_TX2GENE,
        stringtie_bridge=STRINGTIE_BRIDGE,
        gene_namespace=GENE_NAMESPACE
    output:
        TX2GENE_MASTER
    conda:
        "../../envs/qc.yaml"
    log:
        "logs/reference_namespace/build_tx2gene_master.log"
    shell:
        """
        python workflow/scripts/build_tx2gene_master.py \
            --reference-tx2gene {input.reference_tx2gene} \
            --stringtie-bridge {input.stringtie_bridge} \
            --gene-namespace {input.gene_namespace} \
            --output {output} > {log} 2>&1
        """


rule build_import_manifests:
    """Build manifest files for tximport-based quantifiers."""
    input:
        samples=config["samples"],
        salmon=expand("results/05.quantification/native/salmon/per_sample/{sample}/quant.sf", sample=SAMPLES),
        kallisto=expand("results/05.quantification/native/kallisto/per_sample/{sample}/abundance.tsv", sample=SAMPLES),
        stringtie=expand("results/05.quantification/native/stringtie/per_sample/{sample}/final/t_data.ctab", sample=SAMPLES)
    output:
        salmon=SALMON_MANIFEST,
        kallisto=KALLISTO_MANIFEST,
        stringtie=STRINGTIE_MANIFEST
    conda:
        "../../envs/qc.yaml"
    log:
        "logs/reference_namespace/build_import_manifests.log"
    shell:
        """
        python workflow/scripts/build_import_manifests.py \
            --samples {input.samples} \
            --salmon-output {output.salmon} \
            --kallisto-output {output.kallisto} \
            --stringtie-output {output.stringtie} > {log} 2>&1
        """
