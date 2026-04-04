import pandas as pd
import os
from snakemake.shell import shell

# POSIX shell fix as suggested by the audit
shell.executable("bash")

# Configuration file
configfile: "config/config.yaml"

include: "reference_config.smk"

# Get selected aligner from config
ALIGNER = config.get("aligner", "hisat2")

# Sample information
samples_df = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
SAMPLES = samples_df["sample"].tolist()

wildcard_constraints:
    sample="[^/]+"

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
            "r1": f"results/03.decontam/clean/{wildcards.sample}_R1_clean.fastq.gz",
            "r2": f"results/03.decontam/clean/{wildcards.sample}_R2_clean.fastq.gz",
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
        expand("results/03.decontam/clean/{sample}_R1_clean.fastq.gz", sample=samples),
        expand("results/03.decontam/clean/{sample}_R2_clean.fastq.gz", sample=samples),
        expand("results/03.decontam/stats/{sample}_decision_summary.tsv", sample=samples),
        expand("results/03.decontam/qc/{sample}_{read}_clean_fastqc.html", sample=samples, read=["R1", "R2"]),
        expand("results/03.decontam/qc/{sample}_{read}_clean_fastqc.zip", sample=samples, read=["R1", "R2"]),
        "results/03.decontam/stats/project_decontam_summary.tsv",
        "results/03.decontam/reference/contam_scaffolds_blacklist.tsv",
    ]

def get_decontam_clues_targets(samples):
    if not decontam_enabled():
        return []

    return [
        "results/03.decontam/clues/tables/sample_microbial_burden.tsv",
        "results/03.decontam/clues/tables/priority_targets.tsv",
        "results/03.decontam/clues/plots/microbial_composition_stacked_bar.pdf",
        "results/03.decontam/clues/plots/host_context_overlay.pdf",
    ]

def get_selected_alignment_bam(wildcards):
    if ALIGNER == "hisat2":
        return f"results/04.alignment/hisat2/{wildcards.sample}.bam"
    elif ALIGNER == "star":
        return f"results/04.alignment/star/{wildcards.sample}.bam"
    else:
        raise ValueError(f"Unknown aligner: {ALIGNER}. Choose 'hisat2' or 'star'")

def get_selected_alignment_bai(wildcards):
    if ALIGNER == "hisat2":
        return f"results/04.alignment/hisat2/{wildcards.sample}.bam.bai"
    elif ALIGNER == "star":
        return f"results/04.alignment/star/{wildcards.sample}.bam.bai"
    else:
        raise ValueError(f"Unknown aligner: {ALIGNER}. Choose 'hisat2' or 'star'")

def get_sample_groups():
    return sorted(samples_df["group"].astype(str).drop_duplicates().tolist())

def build_pairwise_contrasts(groups=None):
    resolved_groups = groups if groups is not None else get_sample_groups()
    return [
        f"{left}_vs_{right}"
        for idx, left in enumerate(resolved_groups)
        for right in resolved_groups[idx + 1:]
    ]
