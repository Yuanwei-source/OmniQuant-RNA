#!/usr/bin/env python3
"""Shared helpers for quantification matrix aggregation scripts."""

from pathlib import Path
import pandas as pd
from typing import Dict, List, Callable, Optional, Tuple

from annotation_utils import parse_attributes


def ensure_parent_dir(path: Optional[str]) -> None:
    """Create output parent directory when a path is provided."""
    if not path:
        return
    Path(path).parent.mkdir(parents=True, exist_ok=True)


def load_sample_tables(
    input_dir: str, 
    samples: List[str], 
    relative_filename: str, 
    reader_kwargs: Optional[Dict] = None
) -> Tuple[List[str], Dict[str, pd.DataFrame]]:
    """Load per-sample tables and return (sample_order, tables_by_sample)."""
    kwargs = reader_kwargs or {}
    sample_order = []
    tables = {}

    for sample in samples:
        sample_file = Path(input_dir) / sample / relative_filename
        if not sample_file.exists():
            print(f"Warning: {sample_file} does not exist, skipping {sample}")
            continue

        tables[sample] = pd.read_csv(sample_file, **kwargs)
        sample_order.append(sample)

    return sample_order, tables


def build_matrix(
    tables_by_sample: Dict[str, pd.DataFrame], 
    sample_order: List[str], 
    id_columns: List[str], 
    value_extractor: Callable[[pd.DataFrame], pd.Series]
) -> pd.DataFrame:
    """Build a wide matrix using vectorized O(1) concat operations."""
    if not sample_order:
        raise ValueError("No valid sample tables found")

    series_list = []
    
    for sample in sample_order:
        table = tables_by_sample[sample]
        values = value_extractor(table)
        
        # Ensure we have the same length
        if len(values) != len(table):
            raise ValueError(f"Value length mismatch for sample {sample}")
            
        # Set index to id_columns so we can align on them during concat
        if len(id_columns) == 1:
            index = table[id_columns[0]]
        else:
            index = pd.MultiIndex.from_frame(table[id_columns])
            
        if index.duplicated().any():
            dup_count = index.duplicated().sum()
            raise ValueError(f"Sample {sample} has {dup_count} duplicated feature keys")
            
        s = pd.Series(values.to_numpy(), index=index, name=sample)
        series_list.append(s)

    # Perform a single outer join across all samples via concat
    merged = pd.concat(series_list, axis=1, join="outer")
    
    # Restore the ID columns as regular columns
    merged = merged.reset_index()
    
    return merged


def parse_gtf_tx2gene(gtf_file: str) -> Dict[str, str]:
    """Extract transcript_id -> gene_id mapping from a GTF file."""
    tx2gene = {}
    with open(gtf_file, "r") as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue

            attrs = parse_attributes(fields[8])
            transcript_id = attrs.get("transcript_id") or attrs.get("ID")
            gene_id = attrs.get("gene_id") or attrs.get("Parent")
            if transcript_id and gene_id:
                tx2gene[transcript_id] = gene_id

    return tx2gene


def map_transcript_to_gene(transcript_ids: pd.Series, tx2gene: Dict[str, str]) -> pd.Series:
    """Map transcript IDs to genes with version-stripping fallback."""
    mapped = transcript_ids.map(tx2gene)
    missing = mapped.isna()
    if missing.any():
        stripped = transcript_ids[missing].astype(str).str.replace(r"\.[0-9]+$", "", regex=True)
        mapped.loc[missing] = stripped.map(tx2gene)
    return mapped
