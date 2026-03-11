#!/usr/bin/env python3
"""Build the master tx2gene mapping table used by all tximport-based quantifiers."""

import argparse
import csv
from pathlib import Path

import pandas as pd


def normalize_bool(series, default=False, length=0):
    if series is None:
        return pd.Series(["TRUE" if default else "FALSE"] * length)
    default_value = "TRUE" if default else "FALSE"
    return series.fillna(default_value).astype(str).str.upper().replace({"T": "TRUE", "F": "FALSE"})


def main():
    parser = argparse.ArgumentParser(description="Build master tx2gene table")
    parser.add_argument("--reference-tx2gene", required=True, help="Reference tx2gene TSV")
    parser.add_argument("--stringtie-bridge", required=True, help="StringTie bridge TSV")
    parser.add_argument("--gene-namespace", required=True, help="Gene namespace TSV")
    parser.add_argument("--output", required=True, help="Output TSV path")
    args = parser.parse_args()

    ref_df = pd.read_csv(args.reference_tx2gene, sep="\t")
    namespace_df = pd.read_csv(args.gene_namespace, sep="\t")
    namespace = pd.DataFrame(namespace_df[["gene_id", "gene_name", "allow_consensus_main"]]).copy()
    namespace.columns = ["gene_id_resolved", "gene_name", "allow_consensus_main"]

    ref_master = ref_df.rename(columns={"gene_id": "gene_id_resolved"}).copy()
    ref_master["quantifier"] = "reference"
    ref_master["gene_id_original"] = ref_master["gene_id_resolved"]
    ref_master["is_novel_only"] = "FALSE"
    ref_master["is_ambiguous"] = "FALSE"
    ref_master["resolution_status"] = "direct_reference"
    ref_master["exclusion_reason"] = ""
    ref_master["mapping_source"] = ref_master.get("mapping_source", "reference_gtf")
    ref_master["allow_consensus_main"] = normalize_bool(ref_master.get("allow_consensus_main"), default=True, length=len(ref_master))
    ref_master["is_reference_gene"] = normalize_bool(ref_master.get("is_reference_gene"), default=True, length=len(ref_master))
    ref_master["gene_name_filled"] = ref_master.get("gene_name_filled", pd.Series([1] * len(ref_master))).fillna(1).astype(int)

    bridge_df = pd.read_csv(args.stringtie_bridge, sep="\t")
    stringtie_master = bridge_df.copy()
    stringtie_master["quantifier"] = "stringtie"
    stringtie_master["gene_id_original"] = stringtie_master["stringtie_gene_id"]
    stringtie_master["gene_name_filled"] = stringtie_master.get("gene_name_filled", pd.Series([1] * len(stringtie_master))).fillna(1).astype(int)
    stringtie_master["allow_consensus_main"] = normalize_bool(stringtie_master.get("allow_consensus_main"), default=False, length=len(stringtie_master))
    stringtie_master["is_reference_gene"] = normalize_bool(stringtie_master.get("is_reference_gene"), default=False, length=len(stringtie_master))
    stringtie_master["is_novel_only"] = normalize_bool(stringtie_master.get("is_novel_only"), default=False, length=len(stringtie_master))
    stringtie_master["is_ambiguous"] = normalize_bool(stringtie_master.get("is_ambiguous"), default=False, length=len(stringtie_master))
    stringtie_master = stringtie_master.drop(columns=[col for col in ["stringtie_gene_id", "reference_gene_id"] if col in stringtie_master.columns])

    stringtie_master = stringtie_master.merge(namespace, on="gene_id_resolved", how="left", suffixes=("", "_namespace"))
    stringtie_master["gene_name"] = stringtie_master["gene_name_namespace"].fillna(stringtie_master["gene_name"]).fillna(stringtie_master["gene_id_original"])
    if "allow_consensus_main_namespace" in stringtie_master.columns:
        stringtie_master["allow_consensus_main"] = stringtie_master["allow_consensus_main"].where(
            stringtie_master["allow_consensus_main"] == "TRUE",
            stringtie_master["allow_consensus_main_namespace"].fillna("FALSE")
        )
    stringtie_master = stringtie_master.drop(columns=[col for col in ["gene_name_namespace", "allow_consensus_main_namespace"] if col in stringtie_master.columns])

    master = pd.concat(
        [
            ref_master[[
                "quantifier",
                "transcript_id",
                "gene_id_resolved",
                "gene_id_original",
                "gene_name",
                "gene_name_filled",
                "mapping_source",
                "is_reference_gene",
                "is_novel_only",
                "is_ambiguous",
                "allow_consensus_main",
                "resolution_status",
                "exclusion_reason",
            ]],
            stringtie_master[[
                "quantifier",
                "transcript_id",
                "gene_id_resolved",
                "gene_id_original",
                "gene_name",
                "gene_name_filled",
                "mapping_source",
                "is_reference_gene",
                "is_novel_only",
                "is_ambiguous",
                "allow_consensus_main",
                "resolution_status",
                "exclusion_reason",
            ]],
        ],
        ignore_index=True,
    )

    master = master.drop_duplicates(subset=["quantifier", "transcript_id", "gene_id_resolved", "resolution_status"]).sort_values(["quantifier", "transcript_id"])

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    master.to_csv(output_path, sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)

    print(f"[tx2gene_master] Wrote {len(master)} rows to {output_path}")
    for quantifier, sub_df in master.groupby("quantifier"):
        print(f"[tx2gene_master] {quantifier}: {len(sub_df)} rows")
        for status, count in sub_df["resolution_status"].value_counts().sort_index().items():
            print(f"  - {status}: {count}")


if __name__ == "__main__":
    main()
