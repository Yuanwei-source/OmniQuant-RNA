#!/usr/bin/env python3
"""
Aggregate featureCounts results across all samples
"""

import pandas as pd
import argparse
import os
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description='Aggregate featureCounts results')
    parser.add_argument('--input-dir', required=True, help='Directory containing featureCounts results')
    parser.add_argument('--samples', required=True, nargs='+', help='List of sample names')
    parser.add_argument('--output-counts', required=True, help='Output counts matrix file')
    
    args = parser.parse_args()
    Path(args.output_counts).parent.mkdir(parents=True, exist_ok=True)
    
    # Initialize lists to store data
    counts_data = []
    sample_names = []
    
    # Process each sample
    for sample in args.samples:
        counts_file = Path(args.input_dir) / sample / "counts.txt"
        
        if not counts_file.exists():
            print(f"Warning: {counts_file} does not exist, skipping {sample}")
            continue
            
        # Read featureCounts output (skip first line with program info)
        df = pd.read_csv(counts_file, sep='\t', comment='#', skiprows=1)
        
        if len(counts_data) == 0:
            # First sample - initialize with gene info
            gene_info = df[['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length']].copy()
            counts_data.append(gene_info)
        
        # Extract counts column (last column)
        counts_col = df.iloc[:, -1]  # Last column contains the counts
        counts_col.name = sample
        counts_data.append(counts_col)
        sample_names.append(sample)
    
    if not counts_data:
        print("Error: No valid count files found")
        return
    
    # Combine all data
    combined_df = pd.concat(counts_data, axis=1)
    
    # Save counts matrix
    combined_df.to_csv(args.output_counts, sep='\t', index=False)
    print(f"Created counts matrix with {len(combined_df)} genes and {len(sample_names)} samples")
    print(f"Saved to: {args.output_counts}")

if __name__ == "__main__":
    main()
