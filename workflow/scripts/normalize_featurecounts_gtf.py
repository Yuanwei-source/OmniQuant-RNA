import argparse
import gzip
import os
import re
from collections import Counter
from pathlib import Path
from typing import TextIO, cast
from multiprocessing import Pool

SAFE_ID_RE = re.compile(r"[^A-Za-z0-9._:-]+")

def parse_attributes(raw):
    """Parse GTF/GFF-style attribute field without re.compile."""
    attrs = {}
    for field in str(raw).strip().strip(";").split(";"):
        token = field.strip()
        if not token: continue
        if "=" in token: key, value = token.split("=", 1)
        elif " " in token: key, value = token.split(" ", 1)
        else: continue
        attrs[key.strip()] = value.strip().strip('"')
    return attrs

def open_maybe_gzip(path, mode="rt"):
    handle = gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)
    return cast(TextIO, handle)

def format_attributes(attrs):
    fields = []
    for key, value in attrs.items():
        if value is None or value == "":
            continue
        fields.append('{} "{}"'.format(key, value))
    return "; ".join(fields) + ";"

def strip_gene_prefix(value):
    if not value: return value
    for prefix in ("gene-", "gene:"):
        if value.startswith(prefix): return value[len(prefix):]
    return value

def ensure_gene_id(attrs):
    if attrs.get("gene_id"): return False, "existing"
    t_id = attrs.get("transcript_id") or attrs.get("Parent") or attrs.get("ID")
    g_name = attrs.get("gene_name") or attrs.get("gene")
    if t_id and t_id.startswith(("gene-", "gene:")):
        attrs["gene_id"] = t_id.replace("gene:", "gene-", 1)
        return True, "from_transcript_id"
    if g_name:
        attrs["gene_id"] = g_name if g_name.startswith("gene-") else "gene-{}".format(g_name)
        return True, "from_gene_name"
    if attrs.get("Parent"):
        attrs["gene_id"] = attrs["Parent"]
        return True, "from_parent"
    if attrs.get("ID"):
        attrs["gene_id"] = attrs["ID"]
        return True, "from_id"
    return False, "missing"

def ensure_gene_name(attrs):
    if attrs.get("gene_name"): return False, "existing"
    g_id, t_id = attrs.get("gene_id"), attrs.get("transcript_id") or attrs.get("Parent") or attrs.get("ID")
    if g_id:
        attrs["gene_name"] = strip_gene_prefix(g_id)
        return True, "from_gene_id"
    if t_id and t_id.startswith(("gene-", "gene:")):
        attrs["gene_name"] = strip_gene_prefix(t_id)
        return True, "from_transcript_id"
    return False, "missing"

def ensure_transcript_id(attrs):
    if attrs.get("transcript_id"): return False, "existing"
    if attrs.get("ID"):
        attrs["transcript_id"] = attrs["ID"]
        return True, "from_id"
    if attrs.get("Parent"):
        attrs["transcript_id"] = attrs["Parent"]
        return True, "from_parent"
    return False, "missing"

def ensure_exon_id(attrs, chrom, start, end, strand, line_no):
    if attrs.get("exon_id"): return False, "existing"
    t_id = attrs.get("transcript_id") or attrs.get("Parent")
    g_id = attrs.get("gene_id")
    anchor = t_id or g_id or "line{}".format(line_no)
    safe_anchor = SAFE_ID_RE.sub("_", anchor)
    attrs["exon_id"] = "exon:{}:{}:{}:{}:{}".format(chrom, start, end, strand, safe_anchor)
    return True, "from_coordinates"

def process_chunk(chunk_data):
    start_line_no, lines = chunk_data
    processed_lines = []
    stats = Counter()
    for i, line in enumerate(lines):
        line_no = start_line_no + i
        if line.startswith("#") or not line.strip():
            processed_lines.append(line)
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9:
            processed_lines.append(line)
            continue
        chrom, _, feature, start, end, _, strand, _, attr_str = parts
        attrs = parse_attributes(attr_str)
        original_attrs = dict(attrs)

        stats["lines_total"] += 1
        stats["feature_{}".format(feature)] += 1

        if feature in {"transcript", "exon", "gene"}:
            changed, source = ensure_gene_id(attrs)
            if changed:
                stats["gene_id_added"] += 1
                stats["gene_id_{}".format(source)] += 1
            changed, source = ensure_gene_name(attrs)
            if changed:
                stats["gene_name_added"] += 1
                stats["gene_name_{}".format(source)] += 1

        if feature in {"transcript", "exon"}:
            changed, source = ensure_transcript_id(attrs)
            if changed:
                stats["transcript_id_added"] += 1
                stats["transcript_id_{}".format(source)] += 1

        if feature == "exon":
            changed, source = ensure_exon_id(attrs, chrom, start, end, strand, line_no)
            if changed:
                stats["exon_id_added"] += 1
                stats["exon_id_{}".format(source)] += 1

        if attrs != original_attrs:
            stats["lines_modified"] += 1

        parts[8] = format_attributes(attrs)
        processed_lines.append("\t".join(parts) + "\n")
    return processed_lines, stats

def normalize_gtf(input_path, output_path, summary_path, threads=4):
    overall_stats = Counter()
    output_parent = Path(output_path).parent
    summary_parent = Path(summary_path).parent
    if str(output_parent) != ".": output_parent.mkdir(parents=True, exist_ok=True)
    if str(summary_parent) != ".": summary_parent.mkdir(parents=True, exist_ok=True)

    chunk_size = 50000
    def chunk_generator():
        with open_maybe_gzip(input_path, "rt") as src:
            chunk = []
            start_line_no = 1
            for line_no, line in enumerate(src, 1):
                chunk.append(line)
                if len(chunk) >= chunk_size:
                    yield (start_line_no, chunk)
                    start_line_no = line_no + 1
                    chunk = []
            if chunk:
                yield (start_line_no, chunk)

    with pool_context(threads) as pool, open(output_path, "w") as dst:
        for processed_lines, chunk_stats in pool.imap(process_chunk, chunk_generator()):
            dst.writelines(processed_lines)
            overall_stats.update(chunk_stats)

    with open(summary_path, "w") as summary:
        summary.write("metric\tvalue\n")
        for key in sorted(overall_stats):
            summary.write("{}\t{}\n".format(key, overall_stats[key]))

    print("[featureCounts] Normalized annotation written to: {}".format(output_path))
    print("[featureCounts] Summary written to: {}".format(summary_path))

def pool_context(threads):
    return Pool(processes=threads)

def main():
    parser = argparse.ArgumentParser(description="Normalize GTF attributes for featureCounts.")
    parser.add_argument("--input", required=True, help="Input GTF file")
    parser.add_argument("--output", required=True, help="Output normalized GTF file")
    parser.add_argument("--summary", required=True, help="Output summary TSV file")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads for parallel processing")
    args = parser.parse_args()
    normalize_gtf(args.input, args.output, args.summary, args.threads)

if __name__ == "__main__":
    main()
