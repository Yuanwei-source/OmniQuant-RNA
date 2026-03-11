# Quantification Rules
# Transcript quantification using Kallisto and Salmon

KALLISTO_NATIVE_DIR = "results/04.quantification/native/kallisto/per_sample"
KALLISTO_MATRIX_DIR = "results/04.quantification/matrices/kallisto"

rule kallisto_index:
    """
    Build Kallisto index from transcriptome
    """
    input:
        config["reference"]["transcriptome"]
    output:
        "data/reference/kallisto_index/transcriptome.idx"
    conda:
        "../../envs/quantification.yaml"
    log:
        "logs/kallisto_index.log"
    shell:
        "kallisto index -i {output} {input} 2> {log}"

rule kallisto_quant:
    """
    Quantify transcripts using Kallisto
    """
    input:
        index="data/reference/kallisto_index/transcriptome.idx",
        r1="results/02.trimmed_data/{sample}_R1_trimmed.fastq.gz",
        r2="results/02.trimmed_data/{sample}_R2_trimmed.fastq.gz"
    output:
        abundance=f"{KALLISTO_NATIVE_DIR}" + "/{sample}/abundance.tsv",
        h5=f"{KALLISTO_NATIVE_DIR}" + "/{sample}/abundance.h5",
        run_info=f"{KALLISTO_NATIVE_DIR}" + "/{sample}/run_info.json"
    conda:
        "../../envs/quantification.yaml"
    params:
        outdir=f"{KALLISTO_NATIVE_DIR}" + "/{sample}",
        bootstrap=config.get("kallisto", {}).get("bootstrap", 100),
        extra=config.get("kallisto", {}).get("extra", "")
    log:
        "logs/kallisto/{sample}.log"
    threads: 8
    shell:
        """
        mkdir -p {params.outdir}
        kallisto quant -i {input.index} \
        -o {params.outdir} \
        -b {params.bootstrap} \
        -t {threads} \
        {params.extra} \
        {input.r1} {input.r2} 2> {log}
        """

rule aggregate_kallisto_summary:
    """
    Aggregate Kallisto results across all samples into canonical matrices
    """
    input:
        expand(f"{KALLISTO_NATIVE_DIR}" + "/{sample}/abundance.tsv", sample=SAMPLES)
    output:
        counts=f"{KALLISTO_MATRIX_DIR}/kallisto_transcript_counts_matrix.tsv",
        tpm=f"{KALLISTO_MATRIX_DIR}/kallisto_transcript_tpm_matrix.tsv"
    conda:
        "../../envs/qc.yaml"
    params:
        samples=SAMPLES,
        input_dir=KALLISTO_NATIVE_DIR
    log:
        "logs/kallisto/aggregate_summary.log"
    shell:
        """
        python workflow/scripts/aggregate_kallisto.py \
            --input-dir {params.input_dir} \
            --samples {params.samples} \
            --output-counts {output.counts} \
            --output-tpm {output.tpm} 2> {log}
        """

# Compatibility rule for main workflow
rule kallisto_quantification_results:
    """
    Create symlink for main workflow compatibility
    """
    input:
        f"{KALLISTO_NATIVE_DIR}" + "/{sample}/abundance.tsv"
    output:
        "results/quantification_results/{sample}/abundance.tsv"
    shell:
        "mkdir -p $(dirname {output}) && ln -sf ../../04.quantification/native/kallisto/per_sample/{wildcards.sample}/abundance.tsv {output}"
