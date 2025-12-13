#!/usr/bin/env python3
"""
Aggregate Salmon results across all samples
"""

import pandas as pd
import argparse
import os
import re
from pathlib import Path

def parse_gtf(gtf_file):
    """Extract transcript_id to gene_id mapping from GTF"""
    print(f"Parsing GTF file: {gtf_file}")
    tx2gene = {}
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'): continue
            fields = line.strip().split('\t')
            if len(fields) < 9: continue
            
            attributes = fields[8]
            tx_id_match = re.search(r'transcript_id "([^"]+)"', attributes)
            gene_id_match = re.search(r'gene_id "([^"]+)"', attributes)
            
            if tx_id_match and gene_id_match:
                tx2gene[tx_id_match.group(1)] = gene_id_match.group(1)
    print(f"Found {len(tx2gene)} transcript-to-gene mappings")
    return tx2gene

def main():
    parser = argparse.ArgumentParser(description='Aggregate Salmon quantification results')
    parser.add_argument('--input-dir', required=True, help='Directory containing Salmon results')
    parser.add_argument('--samples', required=True, nargs='+', help='List of sample names')
    parser.add_argument('--output-transcript-counts', required=True, help='Output transcript counts matrix file')
    parser.add_argument('--output-transcript-tpm', required=True, help='Output transcript TPM matrix file')
    parser.add_argument('--output-gene-counts', help='Output gene counts matrix file')
    parser.add_argument('--output-gene-tpm', help='Output gene TPM matrix file')
    parser.add_argument('--gtf', help='GTF file for transcript to gene mapping')
    
    args = parser.parse_args()
    
    # Initialize lists to store data
    counts_data = []
    tpm_data = []
    gene_info = None
    sample_names = []
    
    # Process each sample
    for sample in args.samples:
        quant_file = Path(args.input_dir) / sample / "quant.sf"
        
        if not quant_file.exists():
            print(f"Warning: {quant_file} does not exist, skipping {sample}")
            continue
            
        # Read Salmon quantification output
        df = pd.read_csv(quant_file, sep='\t')
        
        if gene_info is None:
            # First sample - store gene info
            gene_info = df[['Name', 'Length', 'EffectiveLength']].copy()
            gene_info.rename(columns={'Name': 'transcript_id'}, inplace=True)
        
        # Extract counts and TPM
        counts_col = df['NumReads'].copy()
        counts_col.name = sample
        counts_data.append(counts_col)
        
        tpm_col = df['TPM'].copy()
        tpm_col.name = sample
        tpm_data.append(tpm_col)
        
        sample_names.append(sample)
    
    if not counts_data:
        print("Error: No valid quantification files found")
        return
    
    # Combine transcript counts data
    counts_df = pd.concat([gene_info[['transcript_id', 'Length', 'EffectiveLength']]] + counts_data, axis=1)
    counts_df.to_csv(args.output_transcript_counts, sep='\t', index=False)
    
    # Combine transcript TPM data
    tpm_df = pd.concat([gene_info[['transcript_id', 'Length', 'EffectiveLength']]] + tpm_data, axis=1)
    tpm_df.to_csv(args.output_transcript_tpm, sep='\t', index=False)
    
    print(f"Created transcript matrices with {len(gene_info)} transcripts and {len(sample_names)} samples")
    print(f"Transcript Counts saved to: {args.output_transcript_counts}")
    print(f"Transcript TPM saved to: {args.output_transcript_tpm}")

    # Gene level aggregation
    if args.gtf and args.output_gene_counts and args.output_gene_tpm:
        tx2gene = parse_gtf(args.gtf)
        
        # Map transcript IDs to gene IDs
        # Note: gene_info['transcript_id'] contains the transcript IDs from Salmon output
        # We need to handle version numbers if present (e.g., ENST000001.1 vs ENST000001)
        # For now, assume exact match or try stripping version if no match
        
        # Create a mapping series
        tx_ids = counts_df['transcript_id']
        gene_ids = tx_ids.map(tx2gene)
        
        # Check how many mapped
        mapped_count = gene_ids.notna().sum()
        print(f"Mapped {mapped_count} out of {len(tx_ids)} transcripts to genes")
        
        if mapped_count == 0:
            print("Warning: No transcripts mapped to genes. Check ID formats in GTF and Salmon output.")
        else:
            # Add gene_id to dataframes
            counts_df['gene_id'] = gene_ids
            tpm_df['gene_id'] = gene_ids
            
            # Drop unmapped
            counts_df_mapped = counts_df.dropna(subset=['gene_id'])
            tpm_df_mapped = tpm_df.dropna(subset=['gene_id'])
            
            # Aggregate by gene_id
            # For counts: sum
            # For TPM: sum
            # We drop Length and EffectiveLength for gene level as they are transcript specific
            # (or we could calculate weighted average, but sum is standard for counts/TPM)
            
            gene_counts = counts_df_mapped.groupby('gene_id')[sample_names].sum().reset_index()
            gene_tpm = tpm_df_mapped.groupby('gene_id')[sample_names].sum().reset_index()
            
            gene_counts.to_csv(args.output_gene_counts, sep='\t', index=False)
            gene_tpm.to_csv(args.output_gene_tpm, sep='\t', index=False)
            
            print(f"Created gene matrices with {len(gene_counts)} genes")
            print(f"Gene Counts saved to: {args.output_gene_counts}")
            print(f"Gene TPM saved to: {args.output_gene_tpm}")

if __name__ == "__main__":
    main()
