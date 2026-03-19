#!/usr/bin/env python3
"""Combine reference and StringTie mappings into a master tx2gene table using Polars."""

import argparse
from pathlib import Path
import polars as pl


COMMON_COLS = [
    "transcript_id", "gene_id", "gene_name", "gene_name_filled",
    "is_reference_gene", "allow_consensus_main", "mapping_source"
]


def normalize_mapping_schema(df: pl.DataFrame) -> pl.DataFrame:
    return df.select(COMMON_COLS).with_columns([
        pl.col("transcript_id").cast(pl.Utf8),
        pl.col("gene_id").cast(pl.Utf8),
        pl.col("gene_name").cast(pl.Utf8),
        pl.col("gene_name_filled").cast(pl.Int32),
        pl.col("is_reference_gene").cast(pl.Utf8),
        pl.col("allow_consensus_main").cast(pl.Utf8),
        pl.col("mapping_source").cast(pl.Utf8),
    ])


def add_master_metadata(df: pl.DataFrame, quantifier: str, resolution_status: str) -> pl.DataFrame:
    return df.with_columns([
        pl.lit(quantifier).alias("quantifier"),
        pl.col("gene_id").alias("gene_id_original"),
        pl.col("gene_id").alias("gene_id_resolved"),
        pl.lit("FALSE").alias("is_novel_only"),
        pl.lit("FALSE").alias("is_ambiguous"),
        pl.lit(resolution_status).alias("resolution_status"),
    ])


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

    # Combine reference and stringtie bridge mappings with an explicit schema.
    # The bridge can legitimately be empty, so we cannot rely on CSV inference.
    ref_sub = normalize_mapping_schema(ref)
    bridge_sub = normalize_mapping_schema(bridge)

    reference_quant = add_master_metadata(ref_sub, "reference", "reference_exact")

    # StringTie final quantification is performed against the merged GTF, but in
    # this dataset the merged models mostly retain reference transcript IDs.
    # Use the reference mapping as a fallback baseline, then let any explicit
    # bridge rows override matching transcript IDs.
    stringtie_fallback = add_master_metadata(ref_sub, "stringtie", "reference_fallback")
    stringtie_bridge = add_master_metadata(bridge_sub, "stringtie", "stringtie_bridge")
    stringtie_quant = pl.concat([stringtie_bridge, stringtie_fallback], how="vertical").unique(
        subset=["quantifier", "transcript_id"],
        keep="first"
    )

    combined = pl.concat([reference_quant, stringtie_quant], how="vertical")

    # Enrich with namespace information mapping back over resolved details
    # We resolve gene_id explicitly
    combined = combined.join(
        namespace.select(["gene_id", "gene_name"]).unique(subset=["gene_id"]).rename({"gene_name": "ns_gene_name"}),
        on="gene_id",
        how="left"
    )

    combined = combined.with_columns([
        pl.coalesce(["ns_gene_name", "gene_name"]).alias("gene_name")
    ]).drop(["ns_gene_name"])

    out_cols = [
        "quantifier",
        "transcript_id",
        "gene_id",
        "gene_id_original",
        "gene_id_resolved",
        "gene_name",
        "gene_name_filled",
        "is_reference_gene",
        "is_novel_only",
        "is_ambiguous",
        "allow_consensus_main",
        "resolution_status",
        "mapping_source",
    ]
    combined = combined.select(out_cols).sort(["quantifier", "transcript_id"])

    # Ensure output directory exists
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Write output
    combined.write_csv(output_path, separator="\t")
    
    print(f"[tx2gene_master] Built {len(combined)} rows master mapping.")

if __name__ == "__main__":
    main()
