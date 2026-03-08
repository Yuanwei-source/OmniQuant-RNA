# Snakemake workflow for RNA-seq transcriptome quantification
# Author: Yuanwei
# Date: 2025-08-09

import pandas as pd
from snakemake.utils import validate

# Configuration file
configfile: "config/config.yaml"

import os
import subprocess

# Data Sanity Check & Preprocessing: Ensure Reference Fasta has standard IDs (no spaces)
def preprocess_genome_fasta(fasta_path):
    if not fasta_path or not os.path.exists(fasta_path):
        return
    
    # Check if there's any fasta header containing spaces
    try:
        has_spaces = False
        with open(fasta_path, 'r') as f:
            headers_checked = 0
            for line in f:
                if line.startswith('>'):
                    if ' ' in line:
                        has_spaces = True
                        break
                    headers_checked += 1
                    if headers_checked > 10:  # Check up to 10 headers
                        break
        
        if has_spaces:
            print(f"\n[OmniQuant-RNA Info] Reformatting sequence IDs in {fasta_path} (removing spaces and descriptions) to avoid downstream errors with Salmon/FeatureCounts...")
            tmp_path = fasta_path + ".tmp"
            subprocess.run(f"awk '{{print $1}}' {fasta_path} > {tmp_path} && mv {tmp_path} {fasta_path}", shell=True, check=True)
            print("[OmniQuant-RNA Info] Genome FASTA headers cleaned successfully!\n")
    except Exception as e:
        print(f"Warning: Failed to check or format genome fasta: {e}")

preprocess_genome_fasta(config.get("reference", {}).get("genome", ""))


# Get selected aligner from config
ALIGNER = config.get("aligner", "hisat2")

# Sample information
samples_df = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
SAMPLES = samples_df["sample"].tolist()

# Functions to get input files from samples.tsv
def get_r1(wildcards):
    return samples_df.loc[wildcards.sample, "fq1"]

def get_r2(wildcards):
    return samples_df.loc[wildcards.sample, "fq2"]

# Define alignment output paths based on selected aligner
def get_alignment_outputs(samples):
    if ALIGNER == "hisat2":
        return {
            'bam': expand("results/03.alignment/hisat2/{sample}.bam", sample=samples),
            'bai': expand("results/03.alignment/hisat2/{sample}.bam.bai", sample=samples)
        }
    elif ALIGNER == "star":
        return {
            'bam': expand("results/03.alignment/star/{sample}.bam", sample=samples),
            'bai': expand("results/03.alignment/star/{sample}.bam.bai", sample=samples)
        }
    else:
        raise ValueError(f"Unknown aligner: {ALIGNER}. Choose 'hisat2' or 'star'")

ALIGNMENT_OUTPUTS = get_alignment_outputs(SAMPLES)

# Include all modular rules
include: "workflow/rules/annotation_conversion.smk"
include: "workflow/rules/qc.smk"
include: "workflow/rules/alignment.smk"
include: "workflow/rules/quantification_kallisto.smk"
include: "workflow/rules/quantification_salmon.smk"
include: "workflow/rules/quantification_featurecounts.smk"
include: "workflow/rules/quantification_stringtie.smk"
include: "workflow/rules/report.smk"
include: "workflow/rules/dea.smk"

# All target outputs - comprehensive RNA-seq pipeline
rule all:
    input:
        # Annotation format conversion
        "data/reference/annotation_conversion_complete.flag",
        
        # Transcriptome extraction
        config["reference"]["transcriptome"],
        
        # Raw data quality control
        expand("results/01.raw_qc/{sample}_{read}_fastqc.html", sample=SAMPLES, read=["R1", "R2"]),
        expand("results/01.raw_qc/{sample}_{read}_fastqc.zip", sample=SAMPLES, read=["R1", "R2"]),
        
        # Quality trimming with fastp
        expand("results/02.trimmed_data/{sample}_{read}_trimmed.fastq.gz", sample=SAMPLES, read=["R1", "R2"]),
        expand("results/02.trimmed_data/{sample}.fastp.html", sample=SAMPLES),
        
        # Post-trimming quality control
        expand("results/02.trimmed_data/{sample}_{read}_trimmed_fastqc.html", sample=SAMPLES, read=["R1", "R2"]),
        
        # Reference indices
        "data/reference/kallisto_index/transcriptome.idx",
        "data/reference/salmon_index",
        "data/reference/hisat2_index",
        
        # Alignment results
        expand("results/03.alignment/{sample}.bam", sample=SAMPLES),
        expand("results/03.alignment/{sample}.bam.bai", sample=SAMPLES),
        
        # Quantification results - multiple methods
        # Kallisto quantification
        expand("results/04.quantification/kallisto/{sample}/abundance.tsv", sample=SAMPLES),
        expand("results/04.quantification/kallisto/{sample}/abundance.h5", sample=SAMPLES),
        "results/04.quantification/kallisto/all_samples_counts_matrix.txt",
        "results/04.quantification/kallisto/all_samples_tpm_matrix.txt",
        
        # Salmon quantification
        expand("results/04.quantification/salmon/{sample}/quant.sf", sample=SAMPLES),
        "results/04.quantification/salmon/all_samples_transcript_counts_matrix.txt",
        "results/04.quantification/salmon/all_samples_transcript_tpm_matrix.txt",
        "results/04.quantification/salmon/all_samples_gene_counts_matrix.txt",
        "results/04.quantification/salmon/all_samples_gene_tpm_matrix.txt",
        
        # featureCounts quantification
        expand("results/04.quantification/featurecounts/{sample}/counts.txt", sample=SAMPLES),
        "results/04.quantification/featurecounts/all_samples_counts.txt",
        "results/04.quantification/featurecounts/all_samples_counts_matrix.txt",
        
        # StringTie quantification
        expand("results/04.quantification/stringtie/{sample}/final/transcripts.gtf", sample=SAMPLES),
        expand("results/04.quantification/stringtie/{sample}/final/gene_abundances.tab", sample=SAMPLES),
        "results/04.quantification/stringtie/gene_id_mapping.tsv",
        "results/04.quantification/stringtie/all_samples_gene_counts_matrix.txt",
        "results/04.quantification/stringtie/all_samples_transcript_counts_matrix.txt",
        "results/04.quantification/stringtie/all_samples_gene_tpm_matrix.txt",
        "results/04.quantification/stringtie/all_samples_transcript_tpm_matrix.txt",
        "results/04.quantification/stringtie/all_samples_gene_fpkm_matrix.txt",
        "results/04.quantification/stringtie/all_samples_transcript_fpkm_matrix.txt",
        
        # StringTie results with original gene IDs (for downstream analysis)
        "gene_counts_matrix.tsv",
        "transcript_counts_matrix.tsv", 
        "gene_tpm_matrix.tsv",
        "transcript_tpm_matrix.tsv",
        
        # DEA results
        expand("results/05.dea/{quantifier}/integration/PCA_plot.pdf", quantifier=["featurecounts", "stringtie", "salmon", "kallisto"]),
        
        # Reports and summaries
        "results/07.reports/multiqc_report.html"

# Optional rules for different analysis subsets

rule qc_only:
    """Run only quality control steps"""
    input:
        expand("results/01.raw_qc/{sample}_{read}_fastqc.html", sample=SAMPLES, read=["R1", "R2"]),
        expand("results/02.trimmed_data/{sample}_{read}_trimmed.fastq.gz", sample=SAMPLES, read=["R1", "R2"]),
        expand("results/02.trimmed_data/{sample}_{read}_trimmed_fastqc.html", sample=SAMPLES, read=["R1", "R2"])

rule stringtie_only:
    """Run only StringTie quantification with original gene IDs"""
    input:
        expand("results/04.quantification/stringtie/{sample}/final/transcripts.gtf", sample=SAMPLES),
        expand("results/04.quantification/stringtie/{sample}/final/gene_abundances.tab", sample=SAMPLES),
        "results/04.quantification/stringtie/gene_id_mapping.tsv",
        "gene_counts_matrix.tsv",
        "gene_tpm_matrix.tsv"

rule alignment_only:
    """Run quality control and alignment only"""
    input:
        expand("results/02.trimmed_data/{sample}_{read}_trimmed.fastq.gz", sample=SAMPLES, read=["R1", "R2"]),
        expand("results/03.alignment/{sample}.bam", sample=SAMPLES),
        expand("results/03.alignment/{sample}.bam.bai", sample=SAMPLES)

rule quantification_kallisto_only:
    """Run Kallisto quantification only"""
    input:
        "data/reference/kallisto_index/transcriptome.idx",
        expand("results/04.quantification/kallisto/{sample}/abundance.tsv", sample=SAMPLES),
        "results/04.quantification/kallisto/all_samples_counts_matrix.txt",
        "results/04.quantification/kallisto/all_samples_tpm_matrix.txt"

rule quantification_salmon_only:
    """Run Salmon quantification only"""
    input:
        "data/reference/salmon_index",
        expand("results/04.quantification/salmon/{sample}/quant.sf", sample=SAMPLES),
        "results/04.quantification/salmon/all_samples_transcript_counts_matrix.txt",
        "results/04.quantification/salmon/all_samples_transcript_tpm_matrix.txt",
        "results/04.quantification/salmon/all_samples_gene_counts_matrix.txt",
        "results/04.quantification/salmon/all_samples_gene_tpm_matrix.txt"

rule quantification_featurecounts_only:
    """Run featureCounts quantification only"""
    input:
        expand("results/04.quantification/featurecounts/{sample}/counts.txt", sample=SAMPLES),
        "results/04.quantification/featurecounts/all_samples_counts.txt",
        "results/04.quantification/featurecounts/all_samples_counts_matrix.txt"

rule quantification_stringtie_only:
    """Run StringTie quantification only"""
    input:
        expand("results/04.quantification/stringtie/{sample}/final/transcripts.gtf", sample=SAMPLES),
        expand("results/04.quantification/stringtie/{sample}/final/gene_abundances.tab", sample=SAMPLES),
        "results/04.quantification/stringtie/all_samples_gene_counts_matrix.txt",
        "results/04.quantification/stringtie/all_samples_transcript_counts_matrix.txt",
        "results/04.quantification/stringtie/all_samples_gene_tpm_matrix.txt",
        "results/04.quantification/stringtie/all_samples_transcript_tpm_matrix.txt",
        "results/04.quantification/stringtie/all_samples_gene_fpkm_matrix.txt",
        "results/04.quantification/stringtie/all_samples_transcript_fpkm_matrix.txt"
