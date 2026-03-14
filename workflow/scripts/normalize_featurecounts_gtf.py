#!/usr/bin/env python3
import argparse
import gzip
import os
import re
from collections import Counter
from pathlib import Path
from typing import TextIO, cast

from annotation_utils import parse_attributes


SAFE_ID_RE = re.compile(r"[^A-Za-z0-9._:-]+")


def open_maybe_gzip(path, mode="rt"):
    handle = gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)
    return cast(TextIO, handle)


def format_attributes(attrs):
    return "; ".join(
        f'{key} "{value}"'
        for key, value in attrs.items()
        if value is not None and value != ""
    ) + ";"


def strip_gene_prefix(value):
    if not value:
        return value
    for prefix in ("gene-", "gene:"):
        if value.startswith(prefix):
            return value[len(prefix):]
    return value


def ensure_gene_id(attrs):
    if attrs.get("gene_id"):
        return False, "existing"

    transcript_id = attrs.get("transcript_id") or attrs.get("Parent") or attrs.get("ID")
    gene_name = attrs.get("gene_name") or attrs.get("gene")

    if transcript_id and transcript_id.startswith(("gene-", "gene:")):
        attrs["gene_id"] = transcript_id.replace("gene:", "gene-", 1)
        return True, "from_transcript_id"

    if gene_name:
        attrs["gene_id"] = gene_name if gene_name.startswith("gene-") else f"gene-{gene_name}"
        return True, "from_gene_name"

    if attrs.get("Parent"):
        attrs["gene_id"] = attrs["Parent"]
        return True, "from_parent"

    if attrs.get("ID"):
        attrs["gene_id"] = attrs["ID"]
        return True, "from_id"

    return False, "missing"


def ensure_gene_name(attrs):
    if attrs.get("gene_name"):
        return False, "existing"

    gene_id = attrs.get("gene_id")
    transcript_id = attrs.get("transcript_id") or attrs.get("Parent") or attrs.get("ID")

    if gene_id:
        attrs["gene_name"] = strip_gene_prefix(gene_id)
        return True, "from_gene_id"

    if transcript_id and transcript_id.startswith(("gene-", "gene:")):
        attrs["gene_name"] = strip_gene_prefix(transcript_id)
        return True, "from_transcript_id"

    return False, "missing"


def ensure_transcript_id(attrs):
    if attrs.get("transcript_id"):
        return False, "existing"

    if attrs.get("ID"):
        attrs["transcript_id"] = attrs["ID"]
        return True, "from_id"

    if attrs.get("Parent"):
        attrs["transcript_id"] = attrs["Parent"]
        return True, "from_parent"

    return False, "missing"


def ensure_exon_id(attrs, chrom, start, end, strand, line_no):
    if attrs.get("exon_id"):
        return False, "existing"

    transcript_id = attrs.get("transcript_id") or attrs.get("Parent")
    gene_id = attrs.get("gene_id")
    anchor = transcript_id or gene_id or f"line{line_no}"
    safe_anchor = SAFE_ID_RE.sub("_", anchor)
    attrs["exon_id"] = f"exon:{chrom}:{start}:{end}:{strand}:{safe_anchor}"
    return True, "from_coordinates"


def normalize_gtf(input_path, output_path, summary_path):
    stats = Counter()
    output_parent = Path(output_path).parent
    summary_parent = Path(summary_path).parent
    if str(output_parent) != ".":
        output_parent.mkdir(parents=True, exist_ok=True)
    if str(summary_parent) != ".":
        summary_parent.mkdir(parents=True, exist_ok=True)

    with open_maybe_gzip(input_path, "rt") as src, open(output_path, "w") as dst:
        for line_no, line in enumerate(src, 1):
            if line.startswith("#") or not line.strip():
                dst.write(line)
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                dst.write(line)
                continue

            chrom, _, feature, start, end, _, strand, _, attr_str = parts
            attrs = parse_attributes(attr_str)
            original_attrs = dict(attrs)

            stats["lines_total"] += 1
            stats[f"feature_{feature}"] += 1

            if feature in {"transcript", "exon", "gene"}:
                changed, source = ensure_gene_id(attrs)
                if changed:
                    stats["gene_id_added"] += 1
                    stats[f"gene_id_{source}"] += 1

                changed, source = ensure_gene_name(attrs)
                if changed:
                    stats["gene_name_added"] += 1
                    stats[f"gene_name_{source}"] += 1

            if feature in {"transcript", "exon"}:
                changed, source = ensure_transcript_id(attrs)
                if changed:
                    stats["transcript_id_added"] += 1
                    stats[f"transcript_id_{source}"] += 1

            if feature == "exon":
                changed, source = ensure_exon_id(attrs, chrom, start, end, strand, line_no)
                if changed:
                    stats["exon_id_added"] += 1
                    stats[f"exon_id_{source}"] += 1

            if attrs != original_attrs:
                stats["lines_modified"] += 1

            parts[8] = format_attributes(attrs)
            dst.write("\t".join(parts) + "\n")

    with open(summary_path, "w") as summary:
        summary.write("metric\tvalue\n")
        for key in sorted(stats):
            summary.write(f"{key}\t{stats[key]}\n")

    print(f"[featureCounts] Normalized annotation written to: {output_path}")
    print(f"[featureCounts] Summary written to: {summary_path}")
    for key in sorted(stats):
        print(f"{key}\t{stats[key]}")


def main():
    parser = argparse.ArgumentParser(description="Normalize GTF attributes for featureCounts.")
    parser.add_argument("--input", required=True, help="Input GTF file")
    parser.add_argument("--output", required=True, help="Output normalized GTF file")
    parser.add_argument("--summary", required=True, help="Output summary TSV file")
    args = parser.parse_args()

    normalize_gtf(args.input, args.output, args.summary)


if __name__ == "__main__":
    main()
