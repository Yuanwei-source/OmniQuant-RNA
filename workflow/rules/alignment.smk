# Alignment Rules
# Read alignment using HISAT2, STAR, and other aligners

# Rule order to resolve ambiguity
ruleorder: hisat2_align > select_alignment_bam
ruleorder: star_align > select_alignment_bam

# Get selected aligner from config
ALIGNER = config.get("aligner", "hisat2")

# Function to get the appropriate BAM file based on selected aligner
def get_aligner_bam(wildcards):
    if ALIGNER == "hisat2":
        return f"results/03.alignment/hisat2/{wildcards.sample}.bam"
    elif ALIGNER == "star":
        return f"results/03.alignment/star/{wildcards.sample}.bam"
    else:
        raise ValueError(f"Unknown aligner: {ALIGNER}. Choose 'hisat2' or 'star'")

def get_aligner_bai(wildcards):
    if ALIGNER == "hisat2":
        return f"results/03.alignment/hisat2/{wildcards.sample}.bam.bai"
    elif ALIGNER == "star":
        return f"results/03.alignment/star/{wildcards.sample}.bam.bai"
    else:
        raise ValueError(f"Unknown aligner: {ALIGNER}. Choose 'hisat2' or 'star'")

rule hisat2_build_index:
    """
    Build HISAT2 index from reference genome
    """
    input:
        genome=config["reference"]["genome"]
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
        hisat2-build -p {threads} {input.genome} \
        {output}/genome 2> {log}
        """

rule hisat2_align:
    """
    Align reads to reference genome using HISAT2
    """
    input:
        r1="results/02.trimmed_data/{sample}_R1_trimmed.fastq.gz",
        r2="results/02.trimmed_data/{sample}_R2_trimmed.fastq.gz",
        index="data/reference/hisat2_index"
    output:
        bam="results/03.alignment/hisat2/{sample}.bam",
        bai="results/03.alignment/hisat2/{sample}.bam.bai"
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
        genome=config["reference"]["genome"],
        gtf=config["reference"]["gtf"]
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
        STAR --runMode genomeGenerate \
        --genomeDir {output} \
        --genomeFastaFiles {input.genome} \
        --sjdbGTFfile {input.gtf} \
        --sjdbOverhang {params.sjdbOverhang} \
        --runThreadN {threads} 2> {log}
        """

rule star_align:
    """
    Align reads to reference genome using STAR
    """
    input:
        r1="results/02.trimmed_data/{sample}_R1_trimmed.fastq.gz",
        r2="results/02.trimmed_data/{sample}_R2_trimmed.fastq.gz",
        index="data/reference/star_index"
    output:
        bam="results/03.alignment/star/{sample}.bam",
        bai="results/03.alignment/star/{sample}.bam.bai",
        log_final="results/03.alignment/star/{sample}_Log.final.out"
    conda:
        "../../envs/alignment.yaml"
    params:
        outdir="results/03.alignment/star"
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
        bam=get_aligner_bam,
        bai=get_aligner_bai
    output:
        bam="results/03.alignment/{sample}.bam",
        bai="results/03.alignment/{sample}.bam.bai"
    shell:
        """
        ln -sf $(realpath {input.bam}) {output.bam}
        ln -sf $(realpath {input.bai}) {output.bai}
        """
