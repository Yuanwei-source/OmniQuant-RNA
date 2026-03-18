#!/usr/bin/env python3
"""Build a unified reference transcript-to-gene mapping table from a GTF file using Polars."""

import argparse
from pathlib import Path
import polars as pl
import re

def extract_attribute(attr_str: str, key: str) -> str:
    """Helper for Polars string extraction."""
    # Pattern to match strictly 'key "value";' or 'key=value;' etc
    pattern = rf'{key}\s*(?:=|")([^";\n]+)["\n;]?'
    return pattern

def main():
    parser = argparse.ArgumentParser(description="Build transcript-to-gene mapping from reference GTF")
    parser.add_argument("--gtf", required=True, help="Reference GTF file")
    parser.add_argument("--output", required=True, help="Output TSV path")
    args = parser.parse_args()
    
    source_gtf = Path(args.gtf).name
    
    # 1. Read GTF efficiently using polars read_csv (with tab separator) ignoring comments
    df = pl.read_csv(
        args.gtf,
        separator="\t",
        has_header=False,
        comment_prefix="#",
        truncate_ragged_lines=True,
        ignore_errors=True,
        new_columns=["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
    )
    
    # 2. Extract desired features and attributes via regex
    # Common Transcript/Gene Keys
    tx_pat = extract_attribute("attribute", r"(?:transcript_id|ID)")
    gene_pat = extract_attribute("attribute", r"(?:gene_id|Parent)")
    name_pat = extract_attribute("attribute", r"(?:gene_name|gene)")
    
    res = df.with_columns([
        pl.col("attribute").str.extract(tx_pat, 1).alias("transcript_id"),
        pl.col("attribute").str.extract(gene_pat, 1).alias("gene_id"),
        pl.col("attribute").str.extract(name_pat, 1).alias("gene_name")
    ]).filter(
        pl.col("transcript_id").is_not_null() & pl.col("gene_id").is_not_null()
    )
    
    # Fill missing gene name with gene_id
    res = res.with_columns([
        (pl.col("gene_name").is_null()).cast(pl.Int32).alias("gene_name_filled"),
        pl.col("gene_name").fill_null(pl.col("gene_id"))
    ])
    
    res = res.with_columns([
        pl.lit("TRUE").alias("is_reference_gene"),
        pl.lit("TRUE").alias("allow_consensus_main"),
        pl.lit(source_gtf).alias("mapping_source")
    ])
    
    # 3. Handle duplicates: Keep transcript feature primary. Then check conflicts.
    # Group by transcript_id to figure out conflicting gene_id assignments.
    conflicts = res.group_by("transcript_id").agg([
        pl.col("gene_id").n_unique().alias("num_genes")
    ]).filter(pl.col("num_genes") > 1)
    
    valid_res = res.join(conflicts, on="transcript_id", how="anti")
    
    # We want to drop duplicates for the exact same transcript mappings, preferring actual 'transcript' features if possible
    # We can sort by feature so that 'transcript' comes first, then take first
    transcript_features = ["transcript", "mRNA", "mrna", "lnc_RNA", "rRNA", "tRNA"]
    
    valid_res = valid_res.with_columns([
        pl.col("feature").is_in(transcript_features).cast(pl.Int32).alias("is_primary_feat")
    ]).sort(["transcript_id", "is_primary_feat"], descending=[False, True]).unique(subset=["transcript_id"], keep="first")
    
    out_cols = [
        "transcript_id",
        "gene_id",
        "gene_name",
        "gene_name_filled",
        "is_reference_gene",
        "allow_consensus_main",
        "mapping_source"
    ]
    final_df = valid_res.select(out_cols).sort("transcript_id")
    
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    final_df.write_csv(output_path, separator="\t")
    
    print(f"[reference_tx2gene] Wrote {len(final_df)} mappings to {output_path}")
    print(f"[reference_tx2gene] Excluded {len(conflicts)} transcripts with conflicting gene assignments")

if __name__ == "__main__":
    main()
