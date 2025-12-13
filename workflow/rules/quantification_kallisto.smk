# Quantification Rules
# Transcript quantification using Kallisto and Salmon

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
        abundance="results/04.quantification/kallisto/{sample}/abundance.tsv",
        h5="results/04.quantification/kallisto/{sample}/abundance.h5",
        run_info="results/04.quantification/kallisto/{sample}/run_info.json"
    conda:
        "../../envs/quantification.yaml"
    params:
        outdir="results/04.quantification/kallisto/{sample}",
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
    Aggregate Kallisto results across all samples
    """
    input:
        expand("results/04.quantification/kallisto/{sample}/abundance.tsv", sample=SAMPLES)
    output:
        counts="results/04.quantification/kallisto/all_samples_counts_matrix.txt",
        tpm="results/04.quantification/kallisto/all_samples_tpm_matrix.txt"
    conda:
        "../../envs/qc.yaml"
    params:
        samples=SAMPLES,
        input_dir="results/04.quantification/kallisto"
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
        "results/04.quantification/kallisto/{sample}/abundance.tsv"
    output:
        "results/quantification_results/{sample}/abundance.tsv"
    shell:
        "mkdir -p $(dirname {output}) && ln -sf ../../04.quantification/kallisto/{wildcards.sample}/abundance.tsv {output}"
