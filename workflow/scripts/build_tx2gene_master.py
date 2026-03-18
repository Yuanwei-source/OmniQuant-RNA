#!/usr/bin/env python3
"""Combine reference and StringTie mappings into a master tx2gene table using Polars."""

import argparse
from pathlib import Path
import polars as pl


def main():
    parser = argparse.ArgumentParser(description="Build master tx2gene table")
    parser.add_argument("--reference-tx2gene", required=True, help="Reference tx2gene TSV")
    parser.add_argument("--stringtie-bridge", required=True, help="StringTie bridge TSV")
    parser.add_argument("--gene-namespace", required=True, help="Gene namespace TSV")
    parser.add_argument("--output", required=True, help="Output TSV path")
    args = parser.parse_args()

    ref = pl.read_csv(args.reference_tx2gene, separator="\t")
    bridge = pl.read_csv(args.stringtie_bridge, separator="\t")
    namespace = pl.read_csv(args.gene_namespace, separator="\t")

    # Combine reference and stringtie bridge mappings
    common_cols = [
        "transcript_id", "gene_id", "gene_name", "gene_name_filled",
        "is_reference_gene", "allow_consensus_main", "mapping_source"
    ]
    ref_sub = ref.select(common_cols)
    bridge_sub = bridge.select(common_cols)

    combined = pl.concat([ref_sub, bridge_sub], how="vertical")

    # Enrich with namespace information mapping back over resolved details
    # We resolve gene_id explicitly
    combined = combined.join(
        namespace.select(["gene_id", "gene_name"]).rename({"gene_name": "ns_gene_name"}),
        on="gene_id",
        how="left"
    )

    combined = combined.with_columns([
        pl.col("gene_id").alias("gene_id_original"),
        # The bridge maps stringtie ID to reference, so the literal gene_id is the resolved one if it exists
        pl.col("gene_id").alias("gene_id_resolved"),
        pl.coalesce(["ns_gene_name", "gene_name"]).alias("gene_name")
    ]).drop(["ns_gene_name"])

    # Ensure output directory exists
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Write output
    combined.write_csv(output_path, separator="\t")
    
    print(f"[tx2gene_master] Built {len(combined)} rows master mapping.")

if __name__ == "__main__":
    main()
