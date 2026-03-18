#!/usr/bin/env python3
"""Build a conservative StringTie-to-reference mapping bridge using Polars."""

import argparse
from pathlib import Path
import polars as pl
from build_reference_tx2gene import extract_attribute

def main():
    parser = argparse.ArgumentParser(description="Build StringTie bridge")
    parser.add_argument("--merged-gtf", required=True, help="StringTie merged GTF")
    parser.add_argument("--reference-tx2gene", required=True, help="Reference tx2gene TSV")
    parser.add_argument("--output", required=True, help="Output TSV path")
    args = parser.parse_args()

    ref = pl.read_csv(args.reference_tx2gene, separator="\t")
    
    # Fast regex extraction of StringTie GTF using Polars
    # Only reading transcript features speeds this up further
    df = pl.read_csv(
        args.merged_gtf,
        separator="\t",
        has_header=False,
        comment_prefix="#",
        truncate_ragged_lines=True,
        ignore_errors=True,
        new_columns=["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
    ).filter(pl.col("feature") == "transcript")

    tx_pat = extract_attribute("attribute", "transcript_id")
    gene_pat = extract_attribute("attribute", "gene_id")
    ref_gene_pat = extract_attribute("attribute", "ref_gene_id")
    
    extracted = df.with_columns([
        pl.col("attribute").str.extract(tx_pat, 1).alias("transcript_id"),
        pl.col("attribute").str.extract(gene_pat, 1).alias("gene_id"),
        pl.col("attribute").str.extract(ref_gene_pat, 1).alias("ref_gene_id")
    ]).filter(
        pl.col("transcript_id").is_not_null() & 
        (pl.col("gene_id").is_not_null() | pl.col("ref_gene_id").is_not_null())
    )

    # Use ref_gene_id if available, fallback to gene_id
    extracted = extracted.with_columns([
        pl.coalesce(["ref_gene_id", "gene_id"]).alias("resolved_gene")
    ])
    
    # Filter for novel StringTie transcripts (MSTRG.*)
    bridge_candidates = extracted.filter(pl.col("transcript_id").str.starts_with("MSTRG"))
    
    # Only keep those that successfully mapped to a known reference gene in our tx2gene index
    # (Conservative policy)
    ref_genes = ref.select("gene_id").unique()
    valid_bridge = bridge_candidates.join(ref_genes, left_on="resolved_gene", right_on="gene_id", how="inner")
    
    # Join the original names to get canonical formats
    valid_bridge = valid_bridge.join(
        ref.select(["gene_id", "gene_name", "gene_name_filled"]).unique(subset=["gene_id"]),
        left_on="resolved_gene",
        right_on="gene_id",
        how="left"
    )

    valid_bridge = valid_bridge.with_columns([
        pl.lit("FALSE").alias("is_reference_gene"),
        pl.lit("FALSE").alias("allow_consensus_main"),
        pl.lit("stringtie_merge").alias("mapping_source")
    ]).rename({"resolved_gene": "gene_id"})

    out_cols = [
        "transcript_id",
        "gene_id",
        "gene_name",
        "gene_name_filled",
        "is_reference_gene",
        "allow_consensus_main",
        "mapping_source"
    ]
    
    final_df = valid_bridge.select(out_cols).sort("transcript_id")

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    final_df.write_csv(output_path, separator="\t")
    print(f"[stringtie_bridge] Found {len(final_df)} conservative bridge mappings.")

if __name__ == "__main__":
    main()
