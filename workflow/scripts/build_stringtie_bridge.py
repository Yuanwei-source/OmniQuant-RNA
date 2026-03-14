#!/usr/bin/env python3
"""Build a conservative StringTie transcript-to-reference-gene bridge table."""

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

from annotation_utils import parse_attributes


TRANSCRIPT_FEATURES = {"transcript", "mRNA", "mrna"}


def clean_value(value):
    if value is None:
        return ""
    value = str(value).strip()
    if value in {"", ".", "-", "NA", "None"}:
        return ""
    return value


def main():
    parser = argparse.ArgumentParser(description="Build StringTie bridge table")
    parser.add_argument("--merged-gtf", required=True, help="Merged StringTie GTF")
    parser.add_argument("--reference-tx2gene", required=True, help="Reference tx2gene TSV")
    parser.add_argument("--output", required=True, help="Output TSV path")
    args = parser.parse_args()

    ref_df = pd.read_csv(args.reference_tx2gene, sep="\t")
    ref_tx2gene = dict(zip(ref_df["transcript_id"], ref_df["gene_id"]))
    ref_gene_name = dict(zip(ref_df["gene_id"], ref_df["gene_name"]))

    rows = []
    mstrg_to_refs = defaultdict(set)
    ref_to_mstrg = defaultdict(set)

    with open(args.merged_gtf) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] not in TRANSCRIPT_FEATURES:
                continue
            attrs = parse_attributes(parts[8])
            stringtie_gene_id = clean_value(attrs.get("gene_id"))
            transcript_id = clean_value(attrs.get("transcript_id") or attrs.get("ID"))
            if not stringtie_gene_id or not transcript_id:
                continue

            explicit_ref_gene = clean_value(attrs.get("ref_gene_id"))
            fallback_ref_gene = clean_value(ref_tx2gene.get(transcript_id, ""))
            reference_gene_id = explicit_ref_gene or fallback_ref_gene
            if reference_gene_id:
                mstrg_to_refs[stringtie_gene_id].add(reference_gene_id)
                ref_to_mstrg[reference_gene_id].add(stringtie_gene_id)

            raw_gene_name = clean_value(attrs.get("gene_name") or attrs.get("ref_gene_name"))
            gene_name = raw_gene_name or ref_gene_name.get(reference_gene_id, "") or reference_gene_id or stringtie_gene_id
            rows.append(
                {
                    "transcript_id": transcript_id,
                    "stringtie_gene_id": stringtie_gene_id,
                    "reference_gene_id": reference_gene_id,
                    "gene_name": gene_name,
                    "gene_name_filled": int(not raw_gene_name),
                }
            )

    ambiguous_mstrg = {gene_id for gene_id, refs in mstrg_to_refs.items() if len(refs) > 1}

    output_rows = []
    status_counts = Counter()
    exclusion_counts = Counter()
    for row in rows:
        stringtie_gene_id = row["stringtie_gene_id"]
        reference_gene_id = row["reference_gene_id"]
        gene_name = row["gene_name"]
        gene_name_filled = row["gene_name_filled"]

        is_ambiguous = stringtie_gene_id in ambiguous_mstrg
        is_novel_only = not reference_gene_id and not is_ambiguous
        allow_consensus_main = (not is_ambiguous) and bool(reference_gene_id)
        exclusion_reason = ""
        if is_ambiguous:
            resolution_status = "ambiguous_excluded"
            gene_id_resolved = ""
            exclusion_reason = "one_mstrg_multiple_reference_genes"
            exclusion_counts[exclusion_reason] += 1
        elif reference_gene_id:
            resolution_status = "collapsed_to_reference" if len(ref_to_mstrg[reference_gene_id]) > 1 else "direct_reference"
            gene_id_resolved = reference_gene_id
        else:
            resolution_status = "novel_only"
            gene_id_resolved = ""
            exclusion_reason = "no_reference_gene_mapping"
            exclusion_counts[exclusion_reason] += 1

        status_counts[resolution_status] += 1
        output_rows.append(
            {
                "transcript_id": row["transcript_id"],
                "stringtie_gene_id": stringtie_gene_id,
                "reference_gene_id": reference_gene_id,
                "gene_id_resolved": gene_id_resolved,
                "gene_name": gene_name or gene_id_resolved or stringtie_gene_id,
                "gene_name_filled": gene_name_filled,
                "reference_gene_count_for_mstrg": len(mstrg_to_refs.get(stringtie_gene_id, set())),
                "mstrg_count_for_reference_gene": len(ref_to_mstrg.get(reference_gene_id, set())) if reference_gene_id else 0,
                "mapping_source": "stringtie_merged_gtf",
                "is_reference_gene": "TRUE" if gene_id_resolved else "FALSE",
                "is_novel_only": "TRUE" if is_novel_only else "FALSE",
                "is_ambiguous": "TRUE" if is_ambiguous else "FALSE",
                "allow_consensus_main": "TRUE" if allow_consensus_main else "FALSE",
                "resolution_status": resolution_status,
                "exclusion_reason": exclusion_reason,
            }
        )

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "transcript_id",
        "stringtie_gene_id",
        "reference_gene_id",
        "gene_id_resolved",
        "gene_name",
        "gene_name_filled",
        "reference_gene_count_for_mstrg",
        "mstrg_count_for_reference_gene",
        "mapping_source",
        "is_reference_gene",
        "is_novel_only",
        "is_ambiguous",
        "allow_consensus_main",
        "resolution_status",
        "exclusion_reason",
    ]
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in sorted(output_rows, key=lambda x: (x["stringtie_gene_id"], x["transcript_id"])):
            writer.writerow(row)

    print(f"[stringtie_bridge] Wrote {len(output_rows)} transcript mappings to {output_path}")
    print("[stringtie_bridge] Exclusion Summary:")
    for key, value in sorted(exclusion_counts.items()):
        print(f"  - {key}: {value}")
    print("[stringtie_bridge] Resolution Summary:")
    for key, value in sorted(status_counts.items()):
        print(f"  - {key}: {value}")
    print(f"[stringtie_bridge] Ambiguous MSTRG loci excluded: {len(ambiguous_mstrg)}")


if __name__ == "__main__":
    main()
