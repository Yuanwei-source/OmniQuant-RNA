# Quantification Rules
# Transcript quantification using Kallisto and Salmon

SALMON_NATIVE_DIR = "results/04.quantification/native/salmon/per_sample"
SALMON_MATRIX_DIR = "results/04.quantification/matrices/salmon"

rule salmon_index:
    """
    Build Salmon index from transcriptome
    """
    input:
        transcriptome=config["reference"]["transcriptome"],
        genome=config["reference"]["genome"]
    output:
        directory("data/reference/salmon_index")
    conda:
        "../../envs/quantification.yaml"
    log:
        "logs/salmon_index.log"
    threads: 12
    shell:
        """
        mkdir -p {output}
        
        grep '^>' {input.genome} | cut -d ' ' -f 1 | sed 's/>//' | sort | uniq > {output}/decoys.txt

        cat {input.transcriptome} {input.genome} > {output}/gentrome.fa

        salmon index -t {output}/gentrome.fa -d {output}/decoys.txt -i {output} -p {threads} 2> {log}
        """

rule salmon_quant:
    """
    Quantify transcripts using Salmon
    """
    input:
        index="data/reference/salmon_index",
        r1="results/02.trimmed_data/{sample}_R1_trimmed.fastq.gz",
        r2="results/02.trimmed_data/{sample}_R2_trimmed.fastq.gz"
    output:
        f"{SALMON_NATIVE_DIR}" + "/{sample}/quant.sf"
    conda:
        "../../envs/quantification.yaml"
    params:
        outdir=f"{SALMON_NATIVE_DIR}" + "/{sample}",
        extra=config.get("salmon", {}).get("extra", "")
    log:
        "logs/salmon/{sample}.log"
    threads: 8
    shell:
        """
        mkdir -p {params.outdir}
        salmon quant -i {input.index} -l A \
        -1 {input.r1} -2 {input.r2} \
        -o {params.outdir} -p {threads} \
        {params.extra} 2> {log}
        """

rule aggregate_salmon_summary:
    """
    Aggregate Salmon results across all samples into canonical matrices
    """
    input:
        quant_files=expand(f"{SALMON_NATIVE_DIR}" + "/{sample}/quant.sf", sample=SAMPLES),
        gtf=config["reference"]["gtf"]
    output:
        transcript_counts=f"{SALMON_MATRIX_DIR}/salmon_transcript_counts_matrix.tsv",
        transcript_tpm=f"{SALMON_MATRIX_DIR}/salmon_transcript_tpm_matrix.tsv",
        gene_counts=f"{SALMON_MATRIX_DIR}/salmon_gene_counts_matrix.tsv",
        gene_tpm=f"{SALMON_MATRIX_DIR}/salmon_gene_tpm_matrix.tsv"
    conda:
        "../../envs/qc.yaml"
    params:
        samples=SAMPLES,
        input_dir=SALMON_NATIVE_DIR
    log:
        "logs/salmon/aggregate_summary.log"
    shell:
        """
        python workflow/scripts/aggregate_salmon.py \
            --input-dir {params.input_dir} \
            --samples {params.samples} \
            --output-transcript-counts {output.transcript_counts} \
            --output-transcript-tpm {output.transcript_tpm} \
            --output-gene-counts {output.gene_counts} \
            --output-gene-tpm {output.gene_tpm} \
            --gtf {input.gtf} 2> {log}
        """

# Compatibility rule for main workflow
rule quantification_results_salmon:
    """
    Create symlink for main workflow compatibility (Salmon version)
    """
    input:
        f"{SALMON_NATIVE_DIR}" + "/{sample}/quant.sf"
    output:
        "results/quantification_results/{sample}/quant.sf"
    shell:
        "mkdir -p $(dirname {output}) && ln -sf ../../04.quantification/native/salmon/per_sample/{wildcards.sample}/quant.sf {output}"
