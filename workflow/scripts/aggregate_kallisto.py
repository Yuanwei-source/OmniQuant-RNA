#!/usr/bin/env python3
"""
Aggregate Kallisto results across all samples
"""

import pandas as pd
import argparse
import os
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description='Aggregate Kallisto quantification results')
    parser.add_argument('--input-dir', required=True, help='Directory containing Kallisto results')
    parser.add_argument('--samples', required=True, nargs='+', help='List of sample names')
    parser.add_argument('--output-counts', required=True, help='Output counts matrix file')
    parser.add_argument('--output-tpm', required=True, help='Output TPM matrix file')
    
    args = parser.parse_args()
    Path(args.output_counts).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output_tpm).parent.mkdir(parents=True, exist_ok=True)
    
    # Initialize lists to store data
    counts_data = []
    tpm_data = []
    gene_info = None
    sample_names = []
    
    # Process each sample
    for sample in args.samples:
        abundance_file = Path(args.input_dir) / sample / "abundance.tsv"
        
        if not abundance_file.exists():
            print(f"Warning: {abundance_file} does not exist, skipping {sample}")
            continue
            
        # Read Kallisto abundance output
        df = pd.read_csv(abundance_file, sep='\t')
        
        if gene_info is None:
            # First sample - store gene info
            gene_info = df[['target_id', 'length', 'eff_length']].copy()
        
        # Extract counts and TPM
        counts_col = df['est_counts'].copy()
        counts_col.name = sample
        counts_data.append(counts_col)
        
        tpm_col = df['tpm'].copy()
        tpm_col.name = sample
        tpm_data.append(tpm_col)
        
        sample_names.append(sample)
    
    if not counts_data:
        print("Error: No valid abundance files found")
        return
    
    # Combine counts data
    counts_df = pd.concat([gene_info[['target_id', 'length', 'eff_length']]] + counts_data, axis=1)
    counts_df.to_csv(args.output_counts, sep='\t', index=False)
    
    # Combine TPM data
    tpm_df = pd.concat([gene_info[['target_id', 'length', 'eff_length']]] + tpm_data, axis=1)
    tpm_df.to_csv(args.output_tpm, sep='\t', index=False)
    
    print(f"Created matrices with {len(gene_info)} transcripts and {len(sample_names)} samples")
    print(f"Counts saved to: {args.output_counts}")
    print(f"TPM saved to: {args.output_tpm}")

if __name__ == "__main__":
    main()
