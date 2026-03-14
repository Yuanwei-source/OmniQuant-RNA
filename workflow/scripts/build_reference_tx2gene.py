#!/usr/bin/env python3
"""Build a unified reference transcript-to-gene mapping table from a GTF file."""

import argparse
import csv
from collections import defaultdict
from pathlib import Path

from annotation_utils import parse_attributes


TRANSCRIPT_FEATURES = {"transcript", "mRNA", "mrna", "lnc_RNA", "rRNA", "tRNA"}


def main():
    parser = argparse.ArgumentParser(description="Build transcript-to-gene mapping from reference GTF")
    parser.add_argument("--gtf", required=True, help="Reference GTF file")
    parser.add_argument("--output", required=True, help="Output TSV path")
    args = parser.parse_args()

    transcript_rows = {}
    conflicting = defaultdict(set)
    source_gtf = Path(args.gtf).name

    with open(args.gtf) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            feature = parts[2]
            attrs = parse_attributes(parts[8])
            transcript_id = attrs.get("transcript_id") or attrs.get("ID")
            gene_id = attrs.get("gene_id") or attrs.get("Parent")
            if not transcript_id or not gene_id:
                continue
            if feature not in TRANSCRIPT_FEATURES and transcript_id in transcript_rows:
                continue

            gene_name = attrs.get("gene_name") or attrs.get("gene") or gene_id
            gene_name_filled = int(not attrs.get("gene_name"))
            row = {
                "transcript_id": transcript_id,
                "gene_id": gene_id,
                "gene_name": gene_name,
                "gene_name_filled": gene_name_filled,
                "is_reference_gene": "TRUE",
                "allow_consensus_main": "TRUE",
                "mapping_source": source_gtf,
            }
            if transcript_id in transcript_rows and transcript_rows[transcript_id]["gene_id"] != gene_id:
                conflicting[transcript_id].update({transcript_rows[transcript_id]["gene_id"], gene_id})
                continue
            transcript_rows.setdefault(transcript_id, row)

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "transcript_id",
        "gene_id",
        "gene_name",
        "gene_name_filled",
        "is_reference_gene",
        "allow_consensus_main",
        "mapping_source",
    ]
    written = 0
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for transcript_id in sorted(transcript_rows):
            if transcript_id in conflicting:
                continue
            writer.writerow(transcript_rows[transcript_id])
            written += 1

    print(f"[reference_tx2gene] Wrote {written} mappings to {output_path}")
    print(f"[reference_tx2gene] Excluded {len(conflicting)} transcripts with conflicting gene assignments")


if __name__ == "__main__":
    main()
