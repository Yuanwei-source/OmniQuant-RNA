#!/usr/bin/env python3
"""
Aggregate Salmon results across all samples
"""

import argparse

from aggregate_common import (
    ensure_parent_dir,
    load_sample_tables,
    build_matrix,
    parse_gtf_tx2gene,
    map_transcript_to_gene,
)

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
    ensure_parent_dir(args.output_transcript_counts)
    ensure_parent_dir(args.output_transcript_tpm)
    ensure_parent_dir(args.output_gene_counts)
    ensure_parent_dir(args.output_gene_tpm)

    sample_names, sample_tables = load_sample_tables(
        input_dir=args.input_dir,
        samples=args.samples,
        relative_filename="quant.sf",
        reader_kwargs={"sep": "\t"},
    )

    if not sample_names:
        print("Error: No valid quantification files found")
        return

    counts_df = build_matrix(
        tables_by_sample=sample_tables,
        sample_order=sample_names,
        id_columns=["Name", "Length", "EffectiveLength"],
        value_extractor=lambda df: df["NumReads"],
    ).rename(columns={"Name": "transcript_id"})
    counts_df.to_csv(args.output_transcript_counts, sep='\t', index=False)

    tpm_df = build_matrix(
        tables_by_sample=sample_tables,
        sample_order=sample_names,
        id_columns=["Name", "Length", "EffectiveLength"],
        value_extractor=lambda df: df["TPM"],
    ).rename(columns={"Name": "transcript_id"})
    tpm_df.to_csv(args.output_transcript_tpm, sep='\t', index=False)

    print(f"Created transcript matrices with {len(counts_df)} transcripts and {len(sample_names)} samples")
    print(f"Transcript Counts saved to: {args.output_transcript_counts}")
    print(f"Transcript TPM saved to: {args.output_transcript_tpm}")

    # Gene level aggregation
    if args.gtf and args.output_gene_counts and args.output_gene_tpm:
        print(f"Parsing GTF file: {args.gtf}")
        tx2gene = parse_gtf_tx2gene(args.gtf)
        print(f"Found {len(tx2gene)} transcript-to-gene mappings")
        
        # Map transcript IDs to gene IDs
        # Note: gene_info['transcript_id'] contains the transcript IDs from Salmon output
        # We need to handle version numbers if present (e.g., ENST000001.1 vs ENST000001)
        # For now, assume exact match or try stripping version if no match
        
        # Create a mapping series
        tx_ids = counts_df['transcript_id']
        gene_ids = map_transcript_to_gene(tx_ids, tx2gene)
        
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
