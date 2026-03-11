# Quantification Rules
# Transcript quantification using StringTie

STRINGTIE_NATIVE_DIR = "results/04.quantification/native/stringtie"
STRINGTIE_MATRIX_DIR = "results/04.quantification/matrices/stringtie"
STRINGTIE_AUDIT_DIR = "results/04.quantification/audit/stringtie"

rule stringtie_assemble:
    """
    Assemble transcripts and quantify gene/transcript expression using StringTie
    """
    input:
        bam="results/03.alignment/{sample}.bam",
        gtf=config["reference"]["gtf"]
    output:
        gtf=f"{STRINGTIE_NATIVE_DIR}/per_sample" + "/{sample}/assembly/transcripts.gtf",
        abundance=f"{STRINGTIE_NATIVE_DIR}/per_sample" + "/{sample}/assembly/gene_abundances.tab",
        coverage=f"{STRINGTIE_NATIVE_DIR}/per_sample" + "/{sample}/assembly/t_data.ctab",
        transcript_coverage=f"{STRINGTIE_NATIVE_DIR}/per_sample" + "/{sample}/assembly/i_data.ctab"
    conda:
        "../../envs/stringtie.yaml"
    params:
        outdir=f"{STRINGTIE_NATIVE_DIR}/per_sample" + "/{sample}/assembly",
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
        expand(f"{STRINGTIE_NATIVE_DIR}/per_sample" + "/{sample}/assembly/transcripts.gtf", sample=SAMPLES)
    output:
        f"{STRINGTIE_NATIVE_DIR}/merged/gtf_list.txt"
    shell:
        """
        ls {input} > {output}
        """

rule stringtie_merge:
    """
    Merge all sample GTF files using StringTie merge
    """
    input:
        gtf_list=f"{STRINGTIE_NATIVE_DIR}/merged/gtf_list.txt",
        reference_gtf=config["reference"]["gtf"]
    output:
        merged_gtf=f"{STRINGTIE_NATIVE_DIR}/merged/merged.gtf"
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
        merged_gtf=f"{STRINGTIE_NATIVE_DIR}/merged/merged.gtf"
    output:
        gtf=f"{STRINGTIE_NATIVE_DIR}/per_sample" + "/{sample}/final/transcripts.gtf",
        abundance=f"{STRINGTIE_NATIVE_DIR}/per_sample" + "/{sample}/final/gene_abundances.tab",
        coverage=f"{STRINGTIE_NATIVE_DIR}/per_sample" + "/{sample}/final/t_data.ctab",
        transcript_coverage=f"{STRINGTIE_NATIVE_DIR}/per_sample" + "/{sample}/final/i_data.ctab"
    conda:
        "../../envs/stringtie.yaml"
    params:
        outdir=f"{STRINGTIE_NATIVE_DIR}/per_sample" + "/{sample}/final",
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
        merged_gtf=f"{STRINGTIE_NATIVE_DIR}/merged/merged.gtf"
    output:
        mapping=f"{STRINGTIE_AUDIT_DIR}/gene_id_mapping.tsv"
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
    Aggregate StringTie results across all samples into canonical matrices
    """
    input:
        gtf_files=expand(f"{STRINGTIE_NATIVE_DIR}/per_sample" + "/{sample}/final/transcripts.gtf", sample=SAMPLES),
        abundance_files=expand(f"{STRINGTIE_NATIVE_DIR}/per_sample" + "/{sample}/final/gene_abundances.tab", sample=SAMPLES),
        gene_mapping=f"{STRINGTIE_AUDIT_DIR}/gene_id_mapping.tsv"
    output:
        gene_counts=f"{STRINGTIE_MATRIX_DIR}/stringtie_gene_counts_matrix.tsv",
        transcript_counts=f"{STRINGTIE_MATRIX_DIR}/stringtie_transcript_counts_matrix.tsv",
        gene_tpm=f"{STRINGTIE_MATRIX_DIR}/stringtie_gene_tpm_matrix.tsv",
        transcript_tpm=f"{STRINGTIE_MATRIX_DIR}/stringtie_transcript_tpm_matrix.tsv",
        gene_fpkm=f"{STRINGTIE_MATRIX_DIR}/stringtie_gene_fpkm_matrix.tsv",
        transcript_fpkm=f"{STRINGTIE_MATRIX_DIR}/stringtie_transcript_fpkm_matrix.tsv"
    conda:
        "../../envs/qc.yaml"
    params:
        input_dir=f"{STRINGTIE_NATIVE_DIR}/per_sample",
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
    Generate original-ID audit matrices for StringTie output inspection
    """
    input:
        gtf_files=expand(f"{STRINGTIE_NATIVE_DIR}/per_sample" + "/{sample}/final/transcripts.gtf", sample=SAMPLES),
        abundance_files=expand(f"{STRINGTIE_NATIVE_DIR}/per_sample" + "/{sample}/final/gene_abundances.tab", sample=SAMPLES),
        gene_mapping=f"{STRINGTIE_AUDIT_DIR}/gene_id_mapping.tsv"
    output:
        gene_counts=f"{STRINGTIE_MATRIX_DIR}/stringtie_original_gene_counts_matrix.tsv",
        transcript_counts=f"{STRINGTIE_MATRIX_DIR}/stringtie_original_transcript_counts_matrix.tsv",
        gene_tpm=f"{STRINGTIE_MATRIX_DIR}/stringtie_original_gene_tpm_matrix.tsv",
        transcript_tpm=f"{STRINGTIE_MATRIX_DIR}/stringtie_original_transcript_tpm_matrix.tsv"
    conda:
        "../../envs/qc.yaml"
    params:
        input_dir=f"{STRINGTIE_NATIVE_DIR}/per_sample",
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
        f"{STRINGTIE_NATIVE_DIR}/per_sample" + "/{sample}/final/gene_abundances.tab"
    output:
        "results/quantification_results/{sample}/stringtie_abundances.tab"
    shell:
        "mkdir -p $(dirname {output}) && ln -sf ../../04.quantification/native/stringtie/per_sample/{wildcards.sample}/final/gene_abundances.tab {output}"
