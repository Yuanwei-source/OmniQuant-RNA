#!/usr/bin/env python3

import os
import argparse
import sys
import polars as pl

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Create gene ID mapping table from merged GTF file'
    )
    parser.add_argument('--merged-gtf', required=True,
                        help='Merged GTF file from StringTie')
    parser.add_argument('--output', required=True,
                        help='Output mapping file (TSV format)')
    parser.add_argument('--verbose', action='store_true',
                        help='Enable verbose output')
    return parser.parse_args()

def extract_gene_mapping_polars(merged_gtf_file, verbose=False):
    """Extract gene ID mapping from merged GTF file using polars"""
    
    # Read GTF using polars
    df = pl.read_csv(
        merged_gtf_file,
        separator='\t',
        has_header=False,
        comment_prefix='#',
        new_columns=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'],
        schema_overrides={
            'start': pl.Int64,
            'end': pl.Int64,
            'score': pl.Utf8,
            'strand': pl.Utf8,
            'frame': pl.Utf8,
            'attribute': pl.Utf8,
        },
        null_values='.',
        truncate_ragged_lines=True 
    ).filter(pl.col('feature') == 'transcript')
    
    # Extract gene_id and ref_gene_id from attribute column
    mapping_df = df.select(
        stringtie_id=pl.col('attribute').str.extract(r'gene_id "([^"]+)"', 1),
        original_id=pl.col('attribute').str.extract(r'ref_gene_id "([^"]+)"', 1)
    ).drop_nulls().unique()
    
    return mapping_df

def main():
    args = parse_arguments()

    output_dir = os.path.dirname(args.output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    if args.verbose:
        print(f"Reading merged GTF file: {args.merged_gtf}")
    
    # Extract gene mapping
    mapping_df = extract_gene_mapping_polars(args.merged_gtf, args.verbose)
    
    if args.verbose:
        print(f"Found {mapping_df.height} gene mappings")
    
    # Write mapping to file
    mapping_df.rename({"stringtie_id": "StringTie_ID", "original_id": "Original_ID"}).sort("StringTie_ID").write_csv(args.output, separator='\t')
    
    print(f"Gene mapping written to: {args.output}")
    print(f"Total mappings: {mapping_df.height}")

if __name__ == "__main__":
    main()
