# Quantification Rules
# Transcript quantification using StringTie

STRINGTIE_SNAPSHOT_DIR = "results/tables/raw_matrices/stringtie"
STRINGTIE_QC_SNAPSHOT_DIR = "results/tables/qc_snapshots/stringtie"

rule stringtie_assemble:
    """
    Assemble transcripts and quantify gene/transcript expression using StringTie
    """
    input:
        bam="results/03.alignment/{sample}.bam",
        gtf=config["reference"]["gtf"]
    output:
        gtf="results/04.quantification/stringtie/{sample}/assembly/transcripts.gtf",
        abundance="results/04.quantification/stringtie/{sample}/assembly/gene_abundances.tab",
        coverage="results/04.quantification/stringtie/{sample}/assembly/t_data.ctab",
        transcript_coverage="results/04.quantification/stringtie/{sample}/assembly/i_data.ctab"
    conda:
        "../../envs/stringtie.yaml"
    params:
        outdir="results/04.quantification/stringtie/{sample}/assembly",
        extra=config.get("stringtie", {}).get("extra", "-e -B")
    log:
        "logs/stringtie/{sample}_assembly.log"
    threads: 8
    shell:
        """
        mkdir -p {params.outdir}
        
        stringtie {input.bam} \
        -G {input.gtf} \
        -o {output.gtf} \
        -A {output.abundance} \
        -l {wildcards.sample} \
        -p {threads} \
        {params.extra} 2> {log}
        """

rule stringtie_merge_list:
    """
    Create list of GTF files for StringTie merge
    """
    input:
        expand("results/04.quantification/stringtie/{sample}/assembly/transcripts.gtf", sample=SAMPLES)
    output:
        "results/04.quantification/stringtie/gtf_list.txt"
    shell:
        """
        ls {input} > {output}
        """

rule stringtie_merge:
    """
    Merge all sample GTF files using StringTie merge
    """
    input:
        gtf_list="results/04.quantification/stringtie/gtf_list.txt",
        reference_gtf=config["reference"]["gtf"]
    output:
        merged_gtf="results/04.quantification/stringtie/merged.gtf"
    conda:
        "../../envs/stringtie.yaml"
    params:
        extra=config.get("stringtie", {}).get("merge_extra", "-c 0 -T 1.0")
    log:
        "logs/stringtie/merge.log"
    threads: 4
    shell:
        """
        stringtie --merge \
        -G {input.reference_gtf} \
        -o {output.merged_gtf} \
        -p {threads} \
        {params.extra} \
        {input.gtf_list} 2> {log}
        """

rule stringtie_quantify_final:
    """
    Final quantification using merged GTF
    """
    input:
        bam="results/03.alignment/{sample}.bam",
        merged_gtf="results/04.quantification/stringtie/merged.gtf"
    output:
        gtf="results/04.quantification/stringtie/{sample}/final/transcripts.gtf",
        abundance="results/04.quantification/stringtie/{sample}/final/gene_abundances.tab",
        coverage="results/04.quantification/stringtie/{sample}/final/t_data.ctab",
        transcript_coverage="results/04.quantification/stringtie/{sample}/final/i_data.ctab"
    conda:
        "../../envs/stringtie.yaml"
    params:
        outdir="results/04.quantification/stringtie/{sample}/final",
        extra=config.get("stringtie", {}).get("extra", "-e -B")
    log:
        "logs/stringtie/{sample}_final.log"
    threads: 8
    shell:
        """
        mkdir -p {params.outdir}
        
        stringtie {input.bam} \
        -G {input.merged_gtf} \
        -o {output.gtf} \
        -A {output.abundance} \
        -l {wildcards.sample}_final \
        -p {threads} \
        {params.extra} 2> {log}
        """

# Use final quantification results for aggregation
rule create_gene_mapping:
    """
    Create gene ID mapping table from merged GTF file
    """
    input:
        merged_gtf="results/04.quantification/stringtie/merged.gtf"
    output:
        mapping=f"{STRINGTIE_QC_SNAPSHOT_DIR}/gene_id_mapping.tsv"
    conda:
        "../../envs/qc.yaml"
    log:
        "logs/stringtie/create_gene_mapping.log"
    shell:
        """
        python -u workflow/scripts/create_gene_mapping.py \
            --merged-gtf {input.merged_gtf} \
            --output {output.mapping} \
            --verbose > {log} 2>&1
        """

rule aggregate_stringtie_summary:
    """
    Aggregate StringTie results across all samples with gene ID mapping as snapshot matrices
    """
    input:
        gtf_files=expand("results/04.quantification/stringtie/{sample}/final/transcripts.gtf", sample=SAMPLES),
        abundance_files=expand("results/04.quantification/stringtie/{sample}/final/gene_abundances.tab", sample=SAMPLES),
        gene_mapping=f"{STRINGTIE_QC_SNAPSHOT_DIR}/gene_id_mapping.tsv"
    output:
        gene_counts=f"{STRINGTIE_SNAPSHOT_DIR}/all_samples_gene_counts_matrix.txt",
        transcript_counts=f"{STRINGTIE_SNAPSHOT_DIR}/all_samples_transcript_counts_matrix.txt",
        gene_tpm=f"{STRINGTIE_SNAPSHOT_DIR}/all_samples_gene_tpm_matrix.txt",
        transcript_tpm=f"{STRINGTIE_SNAPSHOT_DIR}/all_samples_transcript_tpm_matrix.txt",
        gene_fpkm=f"{STRINGTIE_SNAPSHOT_DIR}/all_samples_gene_fpkm_matrix.txt",
        transcript_fpkm=f"{STRINGTIE_SNAPSHOT_DIR}/all_samples_transcript_fpkm_matrix.txt"
    conda:
        "../../envs/qc.yaml"
    params:
        input_dir="results/04.quantification/stringtie",
        pattern="*/final/transcripts.gtf",
        read_length=config["read_length"]
    log:
        "logs/stringtie/aggregate_summary.log"
    shell:
        """
        python -u workflow/scripts/aggregate_stringtie.py \
            --input_dir {params.input_dir} \
            --pattern "{params.pattern}" \
            --gene-mapping {input.gene_mapping} \
            --output-gene-counts {output.gene_counts} \
            --output-transcript-counts {output.transcript_counts} \
            --output-gene-tpm {output.gene_tpm} \
            --output-transcript-tpm {output.transcript_tpm} \
            --output-gene-fpkm {output.gene_fpkm} \
            --output-transcript-fpkm {output.transcript_fpkm} \
            --length {params.read_length} \
            --verbose > {log} 2>&1
        """

rule aggregate_stringtie_original:
    """
    Generate snapshot matrices with original gene IDs (for comparison)
    """
    input:
        gtf_files=expand("results/04.quantification/stringtie/{sample}/final/transcripts.gtf", sample=SAMPLES),
        abundance_files=expand("results/04.quantification/stringtie/{sample}/final/gene_abundances.tab", sample=SAMPLES),
        gene_mapping=f"{STRINGTIE_QC_SNAPSHOT_DIR}/gene_id_mapping.tsv"
    output:
        gene_counts=f"{STRINGTIE_SNAPSHOT_DIR}/gene_counts_matrix.tsv",
        transcript_counts=f"{STRINGTIE_SNAPSHOT_DIR}/transcript_counts_matrix.tsv",
        gene_tpm=f"{STRINGTIE_SNAPSHOT_DIR}/gene_tpm_matrix.tsv",
        transcript_tpm=f"{STRINGTIE_SNAPSHOT_DIR}/transcript_tpm_matrix.tsv"
    conda:
        "../../envs/qc.yaml"
    params:
        input_dir="results/04.quantification/stringtie",
        pattern="*/final/transcripts.gtf",
        read_length=config["read_length"]
    log:
        "logs/stringtie/aggregate_original.log"
    shell:
        """
        python -u workflow/scripts/aggregate_stringtie.py \
            --input_dir {params.input_dir} \
            --pattern "{params.pattern}" \
            --gene-mapping {input.gene_mapping} \
            --output-gene-counts {output.gene_counts} \
            --output-transcript-counts {output.transcript_counts} \
            --output-gene-tpm {output.gene_tpm} \
            --output-transcript-tpm {output.transcript_tpm} \
            --length {params.read_length} \
            --verbose > {log} 2>&1
        """

# Compatibility rule for main workflow
rule quantification_results_stringtie:
    """
    Create symlink for main workflow compatibility
    """
    input:
        "results/04.quantification/stringtie/{sample}/final/gene_abundances.tab"
    output:
        "results/quantification_results/{sample}/stringtie_abundances.tab"
    shell:
        "mkdir -p $(dirname {output}) && ln -sf ../../04.quantification/stringtie/{wildcards.sample}/final/gene_abundances.tab {output}"
