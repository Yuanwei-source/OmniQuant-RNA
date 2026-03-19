# Alignment Rules
# Read alignment using HISAT2, STAR, and other aligners

# Rule order to resolve ambiguity
ruleorder: hisat2_align > select_alignment_bam
ruleorder: star_align > select_alignment_bam

rule hisat2_build_index:
    """
    Build HISAT2 index from reference genome
    """
    input:
        genome=REFERENCE_GENOME
    output:
        directory("data/reference/hisat2_index")
    conda:
        "../../envs/alignment.yaml"
    log:
        "logs/hisat2_build_index.log"
    threads: 16
    shell:
        """
        mkdir -p {output}
        tmp_genome=$(mktemp --suffix=.fa)
        trap 'rm -f "$tmp_genome"' EXIT
        awk '/^>/ {{sub(/ .*/, "", $0)}} {{print}}' {input.genome} > "$tmp_genome"
        hisat2-build -p {threads} "$tmp_genome" {output}/genome 2> {log}
        """

rule hisat2_align:
    """
    Align reads to reference genome using HISAT2
    """
    input:
        r1=get_analysis_r1,
        r2=get_analysis_r2,
        index="data/reference/hisat2_index"
    output:
        bam="results/04.alignment/hisat2/{sample}.bam",
        bai="results/04.alignment/hisat2/{sample}.bam.bai"
    conda:
        "../../envs/alignment.yaml"
    params:
        library=config["hisat2"]["library"]
    log:
        "logs/hisat2/{sample}.log"
    threads: 16
    shell:
        """
        (hisat2 -x {input.index}/genome -1 {input.r1} -2 {input.r2} -p {threads} | \
        samtools view -@ {threads} -bS - | \
        samtools sort -@ {threads} -O BAM -o {output.bam} - ) >> {log} 2>&1

        samtools index {output.bam} >> {log} 2>&1
        """

rule star_build_index:
    """
    Build STAR index from reference genome
    """
    input:
        genome=REFERENCE_GENOME,
        gtf=REFERENCE_GTF
    output:
        directory("data/reference/star_index")
    conda:
        "../../envs/alignment.yaml"
    log:
        "logs/star_build_index.log"
    threads: 16
    params:
        sjdbOverhang=config.get("star", {}).get("sjdbOverhang", 149)
    shell:
        """
        mkdir -p {output}
        tmp_genome=$(mktemp --suffix=.fa)
        trap 'rm -f "$tmp_genome"' EXIT
        awk '/^>/ {{sub(/ .*/, "", $0)}} {{print}}' {input.genome} > "$tmp_genome"
        STAR --runMode genomeGenerate \
        --genomeDir {output} \
        --genomeFastaFiles "$tmp_genome" \
        --sjdbGTFfile {input.gtf} \
        --sjdbOverhang {params.sjdbOverhang} \
        --runThreadN {threads} 2> {log}
        """

rule star_align:
    """
    Align reads to reference genome using STAR
    """
    input:
        r1=get_analysis_r1,
        r2=get_analysis_r2,
        index="data/reference/star_index"
    output:
        bam="results/04.alignment/star/{sample}.bam",
        bai="results/04.alignment/star/{sample}.bam.bai",
        log_final="results/04.alignment/star/{sample}_Log.final.out"
    conda:
        "../../envs/alignment.yaml"
    params:
        outdir="results/04.alignment/star"
    log:
        "logs/star/{sample}.log"
    threads: 16
    shell:
        """
        mkdir -p {params.outdir}
        STAR --genomeDir {input.index} \
        --readFilesIn {input.r1} {input.r2} \
        --readFilesCommand zcat \
        --runThreadN {threads} \
        --outFileNamePrefix {params.outdir}/{wildcards.sample}_ \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMstrandField intronMotif \
        --outFilterIntronMotifs RemoveNoncanonical 2> {log}
        
        # Move and rename the output BAM file
        mv {params.outdir}/{wildcards.sample}_Aligned.sortedByCoord.out.bam {output.bam}
        
        # Copy the log file
        cp {params.outdir}/{wildcards.sample}_Log.final.out {output.log_final}
        
        # Index the BAM file
        samtools index {output.bam}
        """

# Select final BAM file based on chosen aligner
rule select_alignment_bam:
    """
    Select BAM file from chosen aligner for downstream analysis
    """
    input:
        bam=get_selected_alignment_bam,
        bai=get_selected_alignment_bai
    output:
        bam="results/04.alignment/{sample}.bam",
        bai="results/04.alignment/{sample}.bam.bai"
    shell:
        """
        ln -sf $(realpath {input.bam}) {output.bam}
        ln -sf $(realpath {input.bai}) {output.bai}
        """
