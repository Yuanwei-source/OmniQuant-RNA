# Snakemake workflow for RNA-seq transcriptome quantification
# Author: Yuanwei
# Date: 2025-08-09

import pandas as pd
from snakemake.utils import validate

# Configuration file
configfile: "config/config.yaml"

import os
import subprocess

import sys

# Auto-detect Reference Files: user can drop arbitrary named .fasta or .gff3 into data/reference
def auto_detect_references():
    ref_dir = "data/reference"
    if not os.path.exists(ref_dir):
        return

    # 1. Detect Genome FASTA
    if config["reference"]["genome"] == "data/reference/genome.fasta":
        fasta_files = [f for f in os.listdir(ref_dir) if f.endswith(('.fa', '.fasta', '.fna')) 
                       and f not in ['genome.fasta', 'transcriptome.fasta'] 
                       and not f.endswith('.tmp') and not f.endswith('.clean')]
        
        target_fasta = os.path.join(ref_dir, "genome.fasta")
        
        if len(fasta_files) == 1:
            source = fasta_files[0]
            source_path = os.path.join(ref_dir, source)
            if os.path.islink(target_fasta):
                os.unlink(target_fasta)
            if not os.path.exists(target_fasta):
                os.symlink(source, target_fasta)
                print(f"\n[OmniQuant-RNA Info] Auto-detected genome FASTA: symlinked '{source}' to 'genome.fasta'")
        elif len(fasta_files) > 1:
            print(f"\n[OmniQuant-RNA Error] Multiple FASTA files found in {ref_dir}: {fasta_files}.")
            print("Please keep only ONE reference genome file, or specify the exact name in config/config.yaml.\n")
            sys.exit(1)

    # 2. Detect Annotation GFF/GTF
    if config["reference"]["gff3"] in ["data/reference/genome.gff3", "data/reference/genome.gtf"]:
        anno_files = [f for f in os.listdir(ref_dir) if f.endswith(('.gff', '.gff3', '.gtf')) 
                  and f not in ['genome.gff3', 'genome.gtf', 'annotation.gtf', 'annotation.gff3']
                  and not f.startswith('genome_converted')
                  and not f.endswith('.featurecounts.gtf')]
                      
        if len(anno_files) == 1:
            source = anno_files[0]
            source_path = os.path.join(ref_dir, source)
            target_name = "genome.gtf" if source.endswith('.gtf') else "genome.gff3"
            target_path = os.path.join(ref_dir, target_name)
            
            if os.path.islink(target_path):
                os.unlink(target_path)
            if not os.path.exists(target_path):
                os.symlink(source, target_path)
                print(f"[OmniQuant-RNA Info] Auto-detected annotation file: symlinked '{source}' to '{target_name}'\n")
            
            # Dynamically update config so the workflow knows the correct starting format
            config["reference"]["gff3"] = target_path
        elif len(anno_files) > 1:
            print(f"\n[OmniQuant-RNA Error] Multiple annotation files found in {ref_dir}: {anno_files}.")
            print("Please keep only ONE annotation file, or specify the exact name in config/config.yaml.\n")
            sys.exit(1)

auto_detect_references()

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
            subprocess.run(f"awk '{{print $1}}' {fasta_path} > {tmp_path} && mv {tmp_path} {fasta_path} && rm -f {fasta_path}.fai", shell=True, check=True)
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

FASTQ_SUFFIXES = [
    ".fastq.gz",
    ".fq.gz",
    ".fastq.bz2",
    ".fq.bz2",
    ".fastq",
    ".fq",
]


def get_fastq_suffix(path):
    lower_path = str(path).lower()
    for suffix in FASTQ_SUFFIXES:
        if lower_path.endswith(suffix):
            return suffix
    raise ValueError(
        f"Unsupported FASTQ extension for '{path}'. Supported suffixes: {', '.join(FASTQ_SUFFIXES)}"
    )


def get_fastqc_alias_path(sample, read_label, source_path):
    suffix = get_fastq_suffix(source_path)
    return os.path.join(".snakemake", "fastqc_aliases", sample, f"{sample}_{read_label}{suffix}")


def decontam_enabled():
    return bool(config.get("decontam", {}).get("enabled", False))


def get_analysis_reads(wildcards):
    if decontam_enabled():
        return {
            "r1": f"results/02.5.decontam/clean/{wildcards.sample}_R1_clean.fastq.gz",
            "r2": f"results/02.5.decontam/clean/{wildcards.sample}_R2_clean.fastq.gz",
        }

    return {
        "r1": f"results/02.trimmed_data/{wildcards.sample}_R1_trimmed.fastq.gz",
        "r2": f"results/02.trimmed_data/{wildcards.sample}_R2_trimmed.fastq.gz",
    }


def get_analysis_r1(wildcards):
    return get_analysis_reads(wildcards)["r1"]


def get_analysis_r2(wildcards):
    return get_analysis_reads(wildcards)["r2"]


def get_decontam_all_targets(samples):
    if not decontam_enabled():
        return []

    return [
        expand("results/02.5.decontam/clean/{sample}_R1_clean.fastq.gz", sample=samples),
        expand("results/02.5.decontam/clean/{sample}_R2_clean.fastq.gz", sample=samples),
        expand("results/02.5.decontam/stats/{sample}_decision_summary.tsv", sample=samples),
        expand("results/02.5.decontam/qc/{sample}_{read}_clean_fastqc.html", sample=samples, read=["R1", "R2"]),
        expand("results/02.5.decontam/qc/{sample}_{read}_clean_fastqc.zip", sample=samples, read=["R1", "R2"]),
        "results/02.5.decontam/stats/project_decontam_summary.tsv",
        "results/02.5.decontam/reference/contam_scaffolds_blacklist.tsv",
    ]

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
include: "workflow/rules/decontam.smk"
include: "workflow/rules/alignment.smk"
include: "workflow/rules/reference_namespace.smk"
include: "workflow/rules/quantification_kallisto.smk"
include: "workflow/rules/quantification_salmon.smk"
include: "workflow/rules/quantification_featurecounts.smk"
include: "workflow/rules/quantification_stringtie.smk"
include: "workflow/rules/report.smk"
include: "workflow/rules/differential_expression.smk"
include: "workflow/rules/consensus_expression.smk"

# All target outputs - comprehensive RNA-seq pipeline
rule all:
    input:
        # Unified namespace and import manifests
        REFERENCE_TX2GENE,
        GENE_NAMESPACE,
        STRINGTIE_BRIDGE,
        TX2GENE_MASTER,
        SALMON_MANIFEST,
        KALLISTO_MANIFEST,
        STRINGTIE_MANIFEST,

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

        # Optional decontamination outputs
        get_decontam_all_targets(SAMPLES),
        
        # Reference indices
        "data/reference/kallisto_index/transcriptome.idx",
        "data/reference/salmon_index",
        "data/reference/hisat2_index",
        
        # Alignment results
        expand("results/03.alignment/{sample}.bam", sample=SAMPLES),
        expand("results/03.alignment/{sample}.bam.bai", sample=SAMPLES),
        
        # Quantification results - multiple methods
        # Kallisto quantification
        expand("results/04.quantification/native/kallisto/per_sample/{sample}/abundance.tsv", sample=SAMPLES),
        expand("results/04.quantification/native/kallisto/per_sample/{sample}/abundance.h5", sample=SAMPLES),
        "results/04.quantification/matrices/kallisto/kallisto_transcript_counts_matrix.tsv",
        "results/04.quantification/matrices/kallisto/kallisto_transcript_tpm_matrix.tsv",
        
        # Salmon quantification
        expand("results/04.quantification/native/salmon/per_sample/{sample}/quant.sf", sample=SAMPLES),
        "results/04.quantification/matrices/salmon/salmon_transcript_counts_matrix.tsv",
        "results/04.quantification/matrices/salmon/salmon_transcript_tpm_matrix.tsv",
        "results/04.quantification/matrices/salmon/salmon_gene_counts_matrix.tsv",
        "results/04.quantification/matrices/salmon/salmon_gene_tpm_matrix.tsv",
        
        # featureCounts quantification
        "data/reference/genome.featurecounts.gtf",
        "data/reference/genome.featurecounts.summary.tsv",
        expand("results/04.quantification/native/featurecounts/per_sample/{sample}/counts.txt", sample=SAMPLES),
        "results/04.quantification/matrices/featurecounts/featurecounts_gene_counts_matrix.tsv",
        
        # StringTie quantification
        expand("results/04.quantification/native/stringtie/per_sample/{sample}/final/transcripts.gtf", sample=SAMPLES),
        expand("results/04.quantification/native/stringtie/per_sample/{sample}/final/gene_abundances.tab", sample=SAMPLES),
        "results/04.quantification/native/stringtie/merged/merged.gtf",
        "results/04.quantification/audit/stringtie/gene_id_mapping.tsv",
        "results/04.quantification/matrices/stringtie/stringtie_gene_counts_matrix.tsv",
        "results/04.quantification/matrices/stringtie/stringtie_transcript_counts_matrix.tsv",
        "results/04.quantification/matrices/stringtie/stringtie_gene_tpm_matrix.tsv",
        "results/04.quantification/matrices/stringtie/stringtie_transcript_tpm_matrix.tsv",
        "results/04.quantification/matrices/stringtie/stringtie_gene_fpkm_matrix.tsv",
        "results/04.quantification/matrices/stringtie/stringtie_transcript_fpkm_matrix.tsv",
        
        # StringTie results with original gene IDs (for downstream analysis)
        "results/04.quantification/matrices/stringtie/stringtie_original_gene_counts_matrix.tsv",
        "results/04.quantification/matrices/stringtie/stringtie_original_transcript_counts_matrix.tsv",
        "results/04.quantification/matrices/stringtie/stringtie_original_gene_tpm_matrix.tsv",
        "results/04.quantification/matrices/stringtie/stringtie_original_transcript_tpm_matrix.tsv",
        
        # DEA results
        expand("results/05.differential_expression/{quantifier}/integration/PCA_plot.pdf", quantifier=["featurecounts", "stringtie", "salmon", "kallisto"]),

        # Consensus DEA results
        expand("results/06.consensus_expression/{contrast}/consensus_results.tsv", contrast=CONSENSUS_CONTRASTS),
        expand("results/06.consensus_expression/{contrast}/consensus_summary.tsv", contrast=CONSENSUS_CONTRASTS),
        expand("results/06.consensus_expression/{contrast}/tier_diagnostics.tsv", contrast=CONSENSUS_CONTRASTS),
        expand("results/06.consensus_expression/{contrast}/significance_membership.tsv", contrast=CONSENSUS_CONTRASTS),
        expand("results/06.consensus_expression/{contrast}/sensitivity_analysis.tsv", contrast=CONSENSUS_CONTRASTS),
        expand("results/06.consensus_expression/{contrast}/logFC_scatter_salmon_vs_featurecounts.pdf", contrast=CONSENSUS_CONTRASTS),
        expand("results/06.consensus_expression/{contrast}/consensus_volcano.pdf", contrast=CONSENSUS_CONTRASTS),
        expand("results/06.consensus_expression/{contrast}/significance_upset.pdf", contrast=CONSENSUS_CONTRASTS),
        
        # Reports and summaries
        "results/07.reports/multiqc_report.html"

# Optional rules for different analysis subsets

rule qc_only:
    """Run only quality control steps"""
    input:
        expand("results/01.raw_qc/{sample}_{read}_fastqc.html", sample=SAMPLES, read=["R1", "R2"]),
        expand("results/02.trimmed_data/{sample}_{read}_trimmed.fastq.gz", sample=SAMPLES, read=["R1", "R2"]),
        expand("results/02.trimmed_data/{sample}_{read}_trimmed_fastqc.html", sample=SAMPLES, read=["R1", "R2"]),
        get_decontam_all_targets(SAMPLES)

rule stringtie_only:
    """Run only StringTie quantification with original gene IDs"""
    input:
        expand("results/04.quantification/native/stringtie/per_sample/{sample}/final/transcripts.gtf", sample=SAMPLES),
        expand("results/04.quantification/native/stringtie/per_sample/{sample}/final/gene_abundances.tab", sample=SAMPLES),
        "results/04.quantification/native/stringtie/merged/merged.gtf",
        "results/04.quantification/audit/stringtie/gene_id_mapping.tsv",
        "results/04.quantification/matrices/stringtie/stringtie_original_gene_counts_matrix.tsv",
        "results/04.quantification/matrices/stringtie/stringtie_original_gene_tpm_matrix.tsv"

rule alignment_only:
    """Run quality control and alignment only"""
    input:
        get_decontam_all_targets(SAMPLES),
        expand("results/03.alignment/{sample}.bam", sample=SAMPLES),
        expand("results/03.alignment/{sample}.bam.bai", sample=SAMPLES)

rule quantification_kallisto_only:
    """Run Kallisto quantification only"""
    input:
        "data/reference/kallisto_index/transcriptome.idx",
        expand("results/04.quantification/native/kallisto/per_sample/{sample}/abundance.tsv", sample=SAMPLES),
        "results/04.quantification/matrices/kallisto/kallisto_transcript_counts_matrix.tsv",
        "results/04.quantification/matrices/kallisto/kallisto_transcript_tpm_matrix.tsv"

rule quantification_salmon_only:
    """Run Salmon quantification only"""
    input:
        "data/reference/salmon_index",
        expand("results/04.quantification/native/salmon/per_sample/{sample}/quant.sf", sample=SAMPLES),
        "results/04.quantification/matrices/salmon/salmon_transcript_counts_matrix.tsv",
        "results/04.quantification/matrices/salmon/salmon_transcript_tpm_matrix.tsv",
        "results/04.quantification/matrices/salmon/salmon_gene_counts_matrix.tsv",
        "results/04.quantification/matrices/salmon/salmon_gene_tpm_matrix.tsv"

rule consensus_only:
    """Run configured consensus DEA contrasts only"""
    input:
        expand("results/06.consensus_expression/{contrast}/consensus_results.tsv", contrast=CONSENSUS_CONTRASTS)

rule quantification_featurecounts_only:
    """Run featureCounts quantification only"""
    input:
        "data/reference/genome.featurecounts.gtf",
        "data/reference/genome.featurecounts.summary.tsv",
        expand("results/04.quantification/native/featurecounts/per_sample/{sample}/counts.txt", sample=SAMPLES),
        "results/04.quantification/matrices/featurecounts/featurecounts_gene_counts_matrix.tsv"

rule quantification_stringtie_only:
    """Run StringTie quantification only"""
    input:
        expand("results/04.quantification/native/stringtie/per_sample/{sample}/final/transcripts.gtf", sample=SAMPLES),
        expand("results/04.quantification/native/stringtie/per_sample/{sample}/final/gene_abundances.tab", sample=SAMPLES),
        "results/04.quantification/native/stringtie/merged/merged.gtf",
        "results/04.quantification/matrices/stringtie/stringtie_gene_counts_matrix.tsv",
        "results/04.quantification/matrices/stringtie/stringtie_transcript_counts_matrix.tsv",
        "results/04.quantification/matrices/stringtie/stringtie_gene_tpm_matrix.tsv",
        "results/04.quantification/matrices/stringtie/stringtie_transcript_tpm_matrix.tsv",
        "results/04.quantification/matrices/stringtie/stringtie_gene_fpkm_matrix.tsv",
        "results/04.quantification/matrices/stringtie/stringtie_transcript_fpkm_matrix.tsv"
