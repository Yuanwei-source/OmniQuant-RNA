import pandas as pd
import os
import subprocess
import sys

# Configuration file
configfile: "config/config.yaml"

# Sample information
samples_df = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
SAMPLES = samples_df["sample"].tolist()

# Functions to get input files from samples.tsv
def get_r1(wildcards):
    return samples_df.loc[wildcards.sample, "fq1"]

def get_r2(wildcards):
    return samples_df.loc[wildcards.sample, "fq2"]

# Auto-detect Reference Files
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
            if os.path.islink(target_fasta): os.unlink(target_fasta)
            if not os.path.exists(target_fasta):
                os.symlink(source, target_fasta)
                print(f"\n[OmniQuant-RNA] Auto-detected: symlinked '{source}' to 'genome.fasta'")

    # 2. Detect Annotation GFF/GTF
    if config["reference"]["gff3"] in ["data/reference/genome.gff3", "data/reference/genome.gtf"]:
        anno_files = [f for f in os.listdir(ref_dir) if f.endswith(('.gff', '.gff3', '.gtf')) 
                      and f not in ['genome.gff3', 'genome.gtf', 'annotation.gtf', 'annotation.gff3']
                      and not f.startswith('genome_converted')]
        if len(anno_files) == 1:
            source = anno_files[0]
            target_name = "genome.gtf" if source.endswith('.gtf') else "genome.gff3"
            target_path = os.path.join(ref_dir, target_name)
            if os.path.islink(target_path): os.unlink(target_path)
            if not os.path.exists(target_path):
                os.symlink(source, target_path)
            config["reference"]["gff3"] = target_path

auto_detect_references()

# Preprocess FASTA headers
def preprocess_genome_fasta(fasta_path):
    if not fasta_path or not os.path.exists(fasta_path): return
    try:
        has_spaces = False
        with open(fasta_path, 'r') as f:
            for i, line in enumerate(f):
                if line.startswith('>'):
                    if ' ' in line: has_spaces = True; break
                if i > 50: break
        if has_spaces:
            tmp_path = fasta_path + ".tmp"
            subprocess.run(f"awk '{{print $1}}' {fasta_path} > {tmp_path} && mv {tmp_path} {fasta_path}", shell=True, check=True)
    except: pass
preprocess_genome_fasta(config.get("reference", {}).get("genome", ""))

# Get selected aligner from config
ALIGNER = config.get("aligner", "hisat2")

include: "rules/annotation_conversion.smk"

