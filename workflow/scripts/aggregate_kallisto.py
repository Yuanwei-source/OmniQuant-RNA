#!/usr/bin/env python3
"""
Aggregate Kallisto results across all samples
"""

import argparse

from aggregate_common import ensure_parent_dir, load_sample_tables, build_matrix

def main():
    parser = argparse.ArgumentParser(description='Aggregate Kallisto quantification results')
    parser.add_argument('--input-dir', required=True, help='Directory containing Kallisto results')
    parser.add_argument('--samples', required=True, nargs='+', help='List of sample names')
    parser.add_argument('--output-counts', required=True, help='Output counts matrix file')
    parser.add_argument('--output-tpm', required=True, help='Output TPM matrix file')
    
    args = parser.parse_args()
    ensure_parent_dir(args.output_counts)
    ensure_parent_dir(args.output_tpm)

    sample_names, sample_tables = load_sample_tables(
        input_dir=args.input_dir,
        samples=args.samples,
        relative_filename="abundance.tsv",
        reader_kwargs={"sep": "\t"},
    )

    if not sample_names:
        print("Error: No valid abundance files found")
        return

    counts_df = build_matrix(
        tables_by_sample=sample_tables,
        sample_order=sample_names,
        id_columns=["target_id", "length", "eff_length"],
        value_extractor=lambda df: df["est_counts"],
    )
    counts_df.to_csv(args.output_counts, sep='\t', index=False)

    tpm_df = build_matrix(
        tables_by_sample=sample_tables,
        sample_order=sample_names,
        id_columns=["target_id", "length", "eff_length"],
        value_extractor=lambda df: df["tpm"],
    )
    tpm_df.to_csv(args.output_tpm, sep='\t', index=False)

    print(f"Created matrices with {len(counts_df)} transcripts and {len(sample_names)} samples")
    print(f"Counts saved to: {args.output_counts}")
    print(f"TPM saved to: {args.output_tpm}")

if __name__ == "__main__":
    main()
