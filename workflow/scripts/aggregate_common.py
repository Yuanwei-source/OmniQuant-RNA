#!/usr/bin/env python3
"""Shared helpers for quantification matrix aggregation scripts."""

from pathlib import Path

import pandas as pd

from annotation_utils import parse_attributes


def ensure_parent_dir(path):
    """Create output parent directory when a path is provided."""
    if not path:
        return
    Path(path).parent.mkdir(parents=True, exist_ok=True)


def load_sample_tables(input_dir, samples, relative_filename, reader_kwargs=None):
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


def build_matrix(tables_by_sample, sample_order, id_columns, value_extractor):
    """Build a wide matrix by explicit ID-key joins across samples."""
    if not sample_order:
        raise ValueError("No valid sample tables found")

    merged = None
    for sample in sample_order:
        table = tables_by_sample[sample]
        values = value_extractor(table)
        if len(values) != len(table):
            raise ValueError(f"Value length mismatch for sample {sample}")

        sample_matrix = table[id_columns].copy()
        sample_matrix[sample] = values.to_numpy()

        duplicated = sample_matrix.duplicated(subset=id_columns, keep=False)
        if duplicated.any():
            duplicated_count = int(duplicated.sum())
            raise ValueError(
                f"Sample {sample} has {duplicated_count} duplicated feature keys for {id_columns}; "
                "cannot safely merge matrices"
            )

        merged = sample_matrix if merged is None else merged.merge(sample_matrix, on=id_columns, how="outer")

    if merged is None:
        raise ValueError("No sample matrices could be constructed")

    return merged


def parse_gtf_tx2gene(gtf_file):
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


def map_transcript_to_gene(transcript_ids, tx2gene):
    """Map transcript IDs to genes with version-stripping fallback."""
    mapped = transcript_ids.map(tx2gene)
    missing = mapped.isna()
    if missing.any():
        stripped = transcript_ids[missing].astype(str).str.replace(r"\.[0-9]+$", "", regex=True)
        mapped.loc[missing] = stripped.map(tx2gene)
    return mapped
