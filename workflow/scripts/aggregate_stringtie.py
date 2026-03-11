#!/usr/bin/env python3

import os
import sys
import re
import csv
import argparse
import glob
from pathlib import Path
from math import ceil
from collections import defaultdict
import warnings


def ensure_parent_dir(path):
    """Create parent directory for an output file if needed."""
    directory = os.path.dirname(path)
    if directory:
        os.makedirs(directory, exist_ok=True)


def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Generate count and expression matrices from StringTie results',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
This script combines functionality from prepDE.py3 and stringtie_expression_matrix.py
to generate both count matrices (for DESeq2) and expression matrices (TPM/FPKM).

Example usage:
  python aggregate_stringtie_fixed.py \\
    --input_dir results/stringtie \\
    --pattern "." \\
    --length 75 \\
    --gene-mapping gene_id_mapping.tsv \\
    --output-gene-counts gene_counts_matrix.tsv \\
    --output-transcript-counts transcript_counts_matrix.tsv \\
    --output-gene-tpm gene_tpm_matrix.tsv \\
    --output-transcript-tpm transcript_tpm_matrix.tsv
        """
    )
    
    # Arguments from prepDE.py3
    parser.add_argument('-i', '--input', '--input_dir', default='.',
                        help="Folder containing all sample sub-directories [default: %(default)s]")
    # Count matrix outputs
    parser.add_argument('--output-gene-counts', default='gene_counts_matrix.tsv',
                        help="Output file for gene count matrix [default: %(default)s]")
    parser.add_argument('--output-transcript-counts', default='transcript_counts_matrix.tsv',
                        help="Output file for transcript count matrix [default: %(default)s]")
    
    # Expression matrix outputs
    parser.add_argument('--output-gene-tpm', default='gene_tpm_matrix.tsv',
                        help="Output file for gene TPM matrix [default: %(default)s]")
    parser.add_argument('--output-transcript-tpm', default='transcript_tpm_matrix.tsv',
                        help="Output file for transcript TPM matrix [default: %(default)s]")
    parser.add_argument('--output-gene-fpkm', default=None,
                        help="Output file for gene FPKM matrix")
    parser.add_argument('--output-transcript-fpkm', default=None,
                        help="Output file for transcript FPKM matrix")
    
    # Other parameters
    parser.add_argument('-l', '--length', default=150, type=int,
                        help="Average read length [default: %(default)s]")
    parser.add_argument('-p', '--pattern', default=".",
                        help="Regular expression that selects sample subdirectories [default: %(default)s]")
    parser.add_argument('-s', '--string', default="MSTRG",
                        help="Prefix used for geneIDs assigned by StringTie [default: %(default)s]")
    parser.add_argument('--gene-mapping', default=None,
                        help="Gene ID mapping file (TSV format) to convert StringTie IDs back to original IDs")
    parser.add_argument('-v', '--verbose', action="store_true",
                        help="Enable verbose processing")
    
    return parser.parse_args()


def load_gene_mapping(mapping_file):
    """Load gene ID mapping from TSV file"""
    gene_mapping = {}
    
    if not mapping_file:
        return gene_mapping
    
    try:
        with open(mapping_file) as f:
            header = f.readline().strip()
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 2:
                    stringtie_id = fields[0]
                    original_id = fields[1]
                    gene_mapping[stringtie_id] = original_id
        
        print(f"Loaded {len(gene_mapping)} gene ID mappings from {mapping_file}")
        
    except FileNotFoundError:
        print(f"Warning: Gene mapping file not found: {mapping_file}")
    except Exception as e:
        print(f"Error loading gene mapping: {e}")
    
    return gene_mapping


def find_samples(input_dir, pattern):
    """Find sample directories and their GTF files"""
    samples = []
    
    if not os.path.isdir(input_dir):
        print(f"Error: directory '{input_dir}' not found!")
        sys.exit(1)
    
    # Use glob to find GTF files based on pattern
    search_pattern = os.path.join(input_dir, pattern.replace("{sample}", "*"))
    gtf_files = glob.glob(search_pattern)
    
    for gtf_path in gtf_files:
        # Extract sample name from path
        rel_path = os.path.relpath(gtf_path, input_dir)
        rel_parts = Path(rel_path).parts
        sample_name = os.path.dirname(rel_path)

        if len(rel_parts) >= 3 and rel_parts[-2] in {"final", "assembly"}:
            sample_name = rel_parts[-3]

        if sample_name:
            samples.append((sample_name, gtf_path))
    
    if len(samples) == 0:
        print(f"Error: no GTF files found under base directory {input_dir}!")
        print(f"Search pattern used: {search_pattern}")
        sys.exit(1)
    
    samples.sort()
    return samples


def get_gene_id(attribute_string, chromosome, transcript_id):
    """Extract gene ID from GTF attributes"""
    RE_GENE_ID = re.compile('gene_id "([^"]+)"')
    RE_GENE_NAME = re.compile('gene_name "([^"]+)"')
    
    r = RE_GENE_ID.search(attribute_string)
    rn = RE_GENE_NAME.search(attribute_string)
    
    if r:
        if rn:
            return r.group(1) + '|' + rn.group(1)
        else:
            return r.group(1)
    return transcript_id


def get_coverage(attribute_string):
    """Extract coverage from GTF attributes"""
    RE_COVERAGE = re.compile(r'cov "([\-\+\d\.]+)"')
    r = RE_COVERAGE.search(attribute_string)
    if r:
        v = float(r.group(1))
        if v < 0.0:
            v = 0.0
        return v
    return 0.0


def get_tpm(attribute_string):
    """Extract TPM from GTF attributes"""
    RE_TPM = re.compile('TPM "([^"]+)"')
    r = RE_TPM.search(attribute_string)
    if r:
        return float(r.group(1))
    return 0.0


def get_fpkm(attribute_string):
    """Extract FPKM from GTF attributes"""
    RE_FPKM = re.compile('FPKM "([^"]+)"')
    r = RE_FPKM.search(attribute_string)
    if r:
        return float(r.group(1))
    return 0.0


def is_transcript(fields):
    """Check if GTF line is a transcript"""
    return len(fields) > 2 and fields[2] == "transcript"


def process_sample_gtf(sample_path, read_length, verbose=False):
    """Process a single sample's GTF file to extract counts and expression values"""
    RE_TRANSCRIPT_ID = re.compile('transcript_id "([^"]+)"')
    
    transcript_counts = {}
    transcript_tpms = {}
    transcript_fpkms = {}
    gene_ids = {}  # transcript_id -> gene_id mapping
    
    if verbose:
        print(f"Processing sample GTF: {sample_path}")
    
    with open(sample_path) as f:
        transcript_len = 0
        coverage = 0
        tpm = 0
        fpkm = 0
        current_transcript_id = None
        current_gene_id = None
        
        for line_no, line in enumerate(f, 1):
            if line.startswith('#'):
                # Check if file was generated with -e option
                if line_no == 1:
                    if '-e' not in line:
                        print(f"Error: sample file {sample_path} was not generated with -e option!")
                        sys.exit(1)
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            if fields[2] == "transcript":
                # Save previous transcript data
                if current_transcript_id and transcript_len > 0:
                    count = int(ceil(coverage * transcript_len / read_length))
                    transcript_counts[current_transcript_id] = count
                    transcript_tpms[current_transcript_id] = tpm
                    transcript_fpkms[current_transcript_id] = fpkm
                
                # Start new transcript
                try:
                    current_transcript_id = RE_TRANSCRIPT_ID.search(fields[8]).group(1)
                    current_gene_id = get_gene_id(fields[8], fields[0], current_transcript_id)
                    coverage = get_coverage(fields[8])
                    tpm = get_tpm(fields[8])
                    fpkm = get_fpkm(fields[8])
                    transcript_len = 0
                    gene_ids[current_transcript_id] = current_gene_id
                except Exception as e:
                    print(f"Problem parsing file {sample_path} at line {line_no}:\n{line}")
                    print(f"Error: {e}")
                    sys.exit(1)
            
            elif fields[2] == "exon" and current_transcript_id:
                # Add exon length
                transcript_len += int(fields[4]) - int(fields[3]) + 1
        
        # Don't forget the last transcript
        if current_transcript_id and transcript_len > 0:
            count = int(ceil(coverage * transcript_len / read_length))
            transcript_counts[current_transcript_id] = count
            transcript_tpms[current_transcript_id] = tpm
            transcript_fpkms[current_transcript_id] = fpkm
    
    return transcript_counts, transcript_tpms, transcript_fpkms, gene_ids


def read_gene_abundances(sample_dir, gene_mapping=None):
    """Read gene abundances from gene_abundances.tab file"""
    gene_abundances_file = os.path.join(sample_dir, "gene_abundances.tab")
    
    if not os.path.exists(gene_abundances_file):
        warnings.warn(f"Gene abundances file not found: {gene_abundances_file}")
        return {}, {}
    
    gene_tpms = {}
    gene_fpkms = {}
    
    with open(gene_abundances_file) as f:
        header = f.readline().strip().split('\t')
        
        # Find TPM and FPKM columns
        tpm_col = fpkm_col = None
        for i, col_name in enumerate(header):
            if 'TPM' in col_name.upper():
                tpm_col = i
            elif 'FPKM' in col_name.upper():
                fpkm_col = i
        
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < len(header):
                continue
            
            gene_id = fields[0]
            
            # Apply gene mapping if available
            if gene_mapping and gene_id in gene_mapping:
                gene_id = gene_mapping[gene_id]
            
            if tpm_col is not None and tpm_col < len(fields):
                try:
                    gene_tpms[gene_id] = float(fields[tpm_col])
                except ValueError:
                    gene_tpms[gene_id] = 0.0
            
            if fpkm_col is not None and fpkm_col < len(fields):
                try:
                    gene_fpkms[gene_id] = float(fields[fpkm_col])
                except ValueError:
                    gene_fpkms[gene_id] = 0.0
    
    return gene_tpms, gene_fpkms


def aggregate_gene_counts(transcript_counts, gene_ids, sample_name, gene_mapping=None):
    """Aggregate transcript counts to gene level"""
    gene_counts = defaultdict(int)
    
    for transcript_id, count in transcript_counts.items():
        gene_id = gene_ids.get(transcript_id, transcript_id)
        
        # Apply gene mapping if available
        if gene_mapping and gene_id in gene_mapping:
            gene_id = gene_mapping[gene_id]
        
        gene_counts[gene_id] += count
    
    return dict(gene_counts)


def write_count_matrix(data_dict, sample_names, output_file, id_type="gene"):
    """Write count matrix in TSV format (compatible with DESeq2)"""
    
    # Get all unique IDs
    all_ids = set()
    for sample_data in data_dict.values():
        all_ids.update(sample_data.keys())
    
    with open(output_file, 'w') as f:
        # Write header
        if id_type == "gene":
            header_prefix = "Gene_ID"
        else:
            header_prefix = "Transcript_ID"
        header = f"{header_prefix}\t" + "\t".join(sample_names) + "\n"
        f.write(header)
        
        # Write data rows
        for feature_id in sorted(all_ids):
            row = [feature_id]
            
            for sample_name in sample_names:
                count = data_dict.get(sample_name, {}).get(feature_id, 0)
                row.append(str(count))
            
            f.write("\t".join(row) + "\n")


def write_expression_matrix(data_dict, sample_names, output_file, id_type="gene"):
    """Write expression matrix in TSV format"""
    
    # Get all unique IDs
    all_ids = set()
    for sample_data in data_dict.values():
        all_ids.update(sample_data.keys())
    
    with open(output_file, 'w') as f:
        # Write header
        if id_type.lower() == "gene":
            header_prefix = "Gene_ID"
        else:
            header_prefix = "Transcript_ID"
        header = f"{header_prefix}\t" + "\t".join(sample_names) + "\n"
        f.write(header)
        
        # Write data rows
        for feature_id in sorted(all_ids):
            row = [feature_id]
            
            for sample_name in sample_names:
                value = data_dict.get(sample_name, {}).get(feature_id, 0.0)
                # Format all numbers with 6 decimal places
                row.append(f"{value:.6f}")
            
            f.write("\t".join(row) + "\n")


def main():
    """Main function"""
    args = parse_arguments()

    output_paths = [
        args.output_gene_counts,
        args.output_transcript_counts,
        args.output_gene_tpm,
        args.output_transcript_tpm,
        args.output_gene_fpkm,
        args.output_transcript_fpkm,
    ]
    for output_path in output_paths:
        if output_path:
            ensure_parent_dir(output_path)
    
    # Load gene mapping if provided
    gene_mapping = load_gene_mapping(args.gene_mapping)
    
    # Find samples
    samples = find_samples(args.input, args.pattern)
    sample_names = [s[0] for s in samples]
    
    if args.verbose:
        print(f"Found {len(samples)} samples:")
        for sample_name, _ in samples:
            print(f"  {sample_name}")
    
    # Initialize data structures
    all_transcript_counts = {}
    all_transcript_tpms = {}
    all_transcript_fpkms = {}
    all_gene_counts = {}
    all_gene_tpms = {}
    all_gene_fpkms = {}
    all_gene_ids = {}
    
    # Process each sample
    for sample_name, gtf_path in samples:
        if args.verbose:
            print(f"\nProcessing sample: {sample_name}")
        
        # Extract data from GTF file
        transcript_counts, transcript_tpms, transcript_fpkms, gene_ids = \
            process_sample_gtf(gtf_path, args.length, args.verbose)
        
        # Store transcript data
        all_transcript_counts[sample_name] = transcript_counts
        all_transcript_tpms[sample_name] = transcript_tpms
        all_transcript_fpkms[sample_name] = transcript_fpkms
        
        # Update gene IDs mapping
        all_gene_ids.update(gene_ids)
        
        # Aggregate gene counts with mapping
        gene_counts = aggregate_gene_counts(transcript_counts, gene_ids, sample_name, gene_mapping)
        all_gene_counts[sample_name] = gene_counts
        
        # Read gene expression values from gene_abundances.tab
        sample_dir = os.path.dirname(gtf_path)
        gene_tpms, gene_fpkms = read_gene_abundances(sample_dir, gene_mapping)
        all_gene_tpms[sample_name] = gene_tpms
        all_gene_fpkms[sample_name] = gene_fpkms
        
        if args.verbose:
            print(f"  Transcripts: {len(transcript_counts)}")
            print(f"  Genes: {len(gene_counts)}")
    
    # Write count matrices (for DESeq2)
    if args.verbose:
        print(f"\nWriting gene count matrix to: {args.output_gene_counts}")
    write_count_matrix(all_gene_counts, sample_names, args.output_gene_counts, "gene")
    
    if args.verbose:
        print(f"Writing transcript count matrix to: {args.output_transcript_counts}")
    write_count_matrix(all_transcript_counts, sample_names, args.output_transcript_counts, "transcript")
    
    # Write TPM expression matrices
    if args.verbose:
        print(f"Writing gene TPM matrix to: {args.output_gene_tpm}")
    write_expression_matrix(all_gene_tpms, sample_names, args.output_gene_tpm, "Gene")
    
    if args.verbose:
        print(f"Writing transcript TPM matrix to: {args.output_transcript_tpm}")
    write_expression_matrix(all_transcript_tpms, sample_names, args.output_transcript_tpm, "Transcript")
    
    # Write FPKM expression matrices if requested
    if args.output_gene_fpkm:
        if args.verbose:
            print(f"Writing gene FPKM matrix to: {args.output_gene_fpkm}")
        write_expression_matrix(all_gene_fpkms, sample_names, args.output_gene_fpkm, "Gene")
    
    if args.output_transcript_fpkm:
        if args.verbose:
            print(f"Writing transcript FPKM matrix to: {args.output_transcript_fpkm}")
        write_expression_matrix(all_transcript_fpkms, sample_names, args.output_transcript_fpkm, "Transcript")
    
    # Print summary
    total_transcripts = len(set().union(*[d.keys() for d in all_transcript_counts.values()]))
    total_genes = len(set().union(*[d.keys() for d in all_gene_counts.values()]))
    
    print(f"\n=== Summary ===")
    print(f"Processed {len(samples)} samples")
    print(f"Total unique transcripts: {total_transcripts}")
    print(f"Total unique genes: {total_genes}")
    print(f"Average read length: {args.length}")
    
    if gene_mapping:
        print(f"Gene mapping applied: {len(gene_mapping)} mappings")
    else:
        print("No gene mapping applied - using StringTie gene IDs")
    
    print(f"\nOutput files:")
    print(f"  Gene count matrix (for DESeq2): {args.output_gene_counts}")
    print(f"  Transcript count matrix: {args.output_transcript_counts}")
    print(f"  Gene TPM matrix: {args.output_gene_tpm}")
    print(f"  Transcript TPM matrix: {args.output_transcript_tpm}")
    
    if args.output_gene_fpkm:
        print(f"  Gene FPKM matrix: {args.output_gene_fpkm}")
    if args.output_transcript_fpkm:
        print(f"  Transcript FPKM matrix: {args.output_transcript_fpkm}")
    
    print("\nDone!")


if __name__ == "__main__":
    main()
