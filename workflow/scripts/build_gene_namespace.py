#!/usr/bin/env python3
"""Build a unified gene namespace table from the reference tx2gene mapping using Polars."""

import argparse
from pathlib import Path
import polars as pl

def main():
    parser = argparse.ArgumentParser(description="Build gene namespace table")
    parser.add_argument("--tx2gene", required=True, help="Reference tx2gene TSV")
    parser.add_argument("--output", required=True, help="Output TSV path")
    args = parser.parse_args()

    df = pl.read_csv(args.tx2gene, separator="\t")
    
    # Fill backwards
    df = df.with_columns([
        pl.when(pl.col("gene_name") == ".").then(pl.col("gene_id")).otherwise(pl.col("gene_name")).fill_null(pl.col("gene_id")).alias("gene_name"),
        pl.col("gene_name_filled").fill_null(1).cast(pl.Int32)
    ])

    namespace = df.group_by("gene_id").agg([
        pl.col("gene_name").first().alias("gene_name"),
        pl.col("gene_name_filled").max().alias("gene_name_filled"),
        pl.col("transcript_id").n_unique().alias("transcript_count")
    ]).sort("gene_id")

    namespace = namespace.with_columns([
        pl.lit("TRUE").alias("is_reference_gene"),
        pl.lit("TRUE").alias("allow_consensus_main"),
        pl.lit("reference_gtf").alias("namespace_source")
    ])

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    namespace.write_csv(args.output, separator="\t")
    
    print(f"[gene_namespace] Wrote {len(namespace)} genes to {output_path}")

if __name__ == "__main__":
    main()
