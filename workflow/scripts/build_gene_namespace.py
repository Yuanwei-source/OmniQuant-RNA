#!/usr/bin/env python3
"""Build a unified gene namespace table from the reference tx2gene mapping."""

import argparse
import csv
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Build gene namespace table")
    parser.add_argument("--tx2gene", required=True, help="Reference tx2gene TSV")
    parser.add_argument("--output", required=True, help="Output TSV path")
    args = parser.parse_args()

    df = pd.read_csv(args.tx2gene, sep="\t")
    df["gene_name"] = df["gene_name"].fillna(df["gene_id"]).replace(".", pd.NA).fillna(df["gene_id"])
    df["gene_name_filled"] = df["gene_name_filled"].fillna(1).astype(int)

    namespace = (
        df.groupby("gene_id", as_index=False)
        .agg(
            gene_name=("gene_name", "first"),
            gene_name_filled=("gene_name_filled", "max"),
            transcript_count=("transcript_id", "nunique"),
        )
        .sort_values("gene_id")
    )
    namespace["is_reference_gene"] = "TRUE"
    namespace["allow_consensus_main"] = "TRUE"
    namespace["namespace_source"] = "reference_gtf"

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    namespace.to_csv(output_path, sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)

    print(f"[gene_namespace] Wrote {len(namespace)} genes to {output_path}")


if __name__ == "__main__":
    main()
