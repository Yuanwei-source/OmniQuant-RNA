#!/usr/bin/env python3
"""
Aggregate featureCounts results across all samples
"""

import argparse

from aggregate_common import ensure_parent_dir, load_sample_tables, build_matrix

def main():
    parser = argparse.ArgumentParser(description='Aggregate featureCounts results')
    parser.add_argument('--input-dir', required=True, help='Directory containing featureCounts results')
    parser.add_argument('--samples', required=True, nargs='+', help='List of sample names')
    parser.add_argument('--output-counts', required=True, help='Output counts matrix file')
    
    args = parser.parse_args()
    ensure_parent_dir(args.output_counts)

    sample_names, sample_tables = load_sample_tables(
        input_dir=args.input_dir,
        samples=args.samples,
        relative_filename="counts.txt",
        reader_kwargs={"sep": "\t", "comment": "#", "skiprows": 1},
    )

    if not sample_names:
        print("Error: No valid count files found")
        return

    combined_df = build_matrix(
        tables_by_sample=sample_tables,
        sample_order=sample_names,
        id_columns=["Geneid", "Chr", "Start", "End", "Strand", "Length"],
        value_extractor=lambda df: df.iloc[:, -1],
    )
    
    # Save counts matrix
    combined_df.to_csv(args.output_counts, sep='\t', index=False)
    print(f"Created counts matrix with {len(combined_df)} genes and {len(sample_names)} samples")
    print(f"Saved to: {args.output_counts}")

if __name__ == "__main__":
    main()
