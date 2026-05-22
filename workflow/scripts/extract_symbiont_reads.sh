#!/usr/bin/env bash
# ============================================================================
# extract_symbiont_reads.sh
#
# Extracts symbiont/microbial read pairs from Kaiju lineage output via seqtk.
# Filters Kaiju lineage by target taxon (e.g., Bacteria;), extracts read IDs
# from column 2, and pulls matching reads from unresolved FASTQ files.
#
# Kaiju lineage format (tab-sep, -p flag):
#   C/U\tread_id\ttaxid\t...\tcellular organisms;Bacteria;Proteobacteria;...
#
# Usage:
#   extract_symbiont_reads.sh \
#     --lineage results/03.decontam/tmp/ND1.kaiju.out.lineage \
#     --r1 results/03.decontam/tmp/ND1_R1_unresolved.fastq.gz \
#     --r2 results/03.decontam/tmp/ND1_R2_unresolved.fastq.gz \
#     --targets "Bacteria;Fungi;Viruses;" \
#     --output-prefix results/03.decontam/clues/extracted/ND1
#
# Output: {prefix}_{target}_R1.fq.gz and {prefix}_{target}_R2.fq.gz
# ============================================================================

set -euo pipefail

DEFAULT_TARGETS="Bacteria;Fungi;Viruses;Wolbachia;"

LINEAGE=""
R1=""
R2=""
TARGETS="$DEFAULT_TARGETS"
OUTPUT_PREFIX=""

usage() {
    cat <<'EOF'
Usage: extract_symbiont_reads.sh [OPTIONS]

Required:
  --lineage FILE       Kaiju lineage output file (tab-separated, -p format)
  --r1 FILE            Unresolved R1 FASTQ (gzipped)
  --r2 FILE            Unresolved R2 FASTQ (gzipped)
  --output-prefix PREFIX  Output file path prefix (e.g., results/.../extracted/ND1)

Optional:
  --targets "LIST"     Semicolon-separated list of target taxa
                       (default: "Bacteria;Fungi;Viruses;Wolbachia;")
  -h, --help           Show this help message

Output:
  {prefix}_{target}_R1.fq.gz
  {prefix}_{target}_R2.fq.gz
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --lineage)       LINEAGE="$2";       shift 2 ;;
        --r1)            R1="$2";            shift 2 ;;
        --r2)            R2="$2";            shift 2 ;;
        --targets)       TARGETS="$2";       shift 2 ;;
        --output-prefix) OUTPUT_PREFIX="$2"; shift 2 ;;
        -h|--help)       usage; exit 0 ;;
        *)
            echo "Error: Unknown option: $1" >&2
            usage >&2; exit 1
            ;;
    esac
done

if [[ -z "$LINEAGE" || -z "$R1" || -z "$R2" || -z "$OUTPUT_PREFIX" ]]; then
    echo "Error: --lineage, --r1, --r2, and --output-prefix are required" >&2
    usage >&2
    exit 1
fi

for f in "$LINEAGE" "$R1" "$R2"; do
    if [[ ! -f "$f" ]]; then
        echo "Error: File not found: $f"
        exit 1
    fi
done

for cmd in seqtk gzip gunzip; do
    if ! command -v "$cmd" &>/dev/null; then
        echo "Error: Required command not found: $cmd"
        exit 1
    fi
done

if command -v pigz &>/dev/null; then
    GZIP_CMD="pigz -c"
else
    GZIP_CMD="gzip -c"
fi

mkdir -p "$(dirname "$OUTPUT_PREFIX")"

TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT

echo "=== extract_symbiont_reads.sh ==="
echo "Lineage:  $LINEAGE"
echo "R1:       $R1"
echo "R2:       $R2"
echo "Targets:  $TARGETS"
echo "Prefix:   $OUTPUT_PREFIX"
echo ""

total_lineage_reads=$(wc -l < "$LINEAGE")
echo "Total reads in lineage file: ${total_lineage_reads}"
echo ""

IFS=';' read -ra TARGET_LIST <<< "$TARGETS"

SUMMARY=""

for target in "${TARGET_LIST[@]}"; do
    [[ -z "$target" ]] && continue

    safe_name=$(echo "$target" | sed 's/[^a-zA-Z0-9._-]/_/g' | sed 's/_*$//')

    echo "--- Processing target: ${target} (safe: ${safe_name}) ---"

    grep -F "${target}" "$LINEAGE" | cut -f2 > "${TMPDIR}/${safe_name}.raw_ids.txt"

    raw_count=$(wc -l < "${TMPDIR}/${safe_name}.raw_ids.txt")

    if [[ "$raw_count" -eq 0 ]]; then
        echo "  No reads found for target '${target}', creating empty output"
        printf '' | $GZIP_CMD > "${OUTPUT_PREFIX}_${safe_name}_R1.fq.gz"
        printf '' | $GZIP_CMD > "${OUTPUT_PREFIX}_${safe_name}_R2.fq.gz"
        SUMMARY="${SUMMARY}  Target '${target}': 0 read pairs\n"
        continue
    fi

    awk '{
        base = $0
        if (base ~ /\/[12]$/) sub(/\/[12]$/, "", base)
        print base
    }' "${TMPDIR}/${safe_name}.raw_ids.txt" | sort -u > "${TMPDIR}/${safe_name}.base_ids.txt"

    awk '{print $0; print $0 "/1"; print $0 "/2"}' \
        "${TMPDIR}/${safe_name}.base_ids.txt" | sort -u > "${TMPDIR}/${safe_name}.ids_for_seqtk.txt"

    unique_pairs=$(wc -l < "${TMPDIR}/${safe_name}.base_ids.txt")
    total_ids=$(wc -l < "${TMPDIR}/${safe_name}.ids_for_seqtk.txt")
    echo "  Raw lineage matches: ${raw_count}, unique read pairs: ${unique_pairs}, seqtk IDs: ${total_ids}"

    seqtk subseq "$R1" "${TMPDIR}/${safe_name}.ids_for_seqtk.txt" | $GZIP_CMD > "${OUTPUT_PREFIX}_${safe_name}_R1.fq.gz"
    seqtk subseq "$R2" "${TMPDIR}/${safe_name}.ids_for_seqtk.txt" | $GZIP_CMD > "${OUTPUT_PREFIX}_${safe_name}_R2.fq.gz"

    extracted_r1_lines=$(gunzip -c "${OUTPUT_PREFIX}_${safe_name}_R1.fq.gz" | wc -l)
    extracted_pairs=$(( extracted_r1_lines / 4 ))

    echo "  Extracted ${extracted_pairs} read pairs for target '${target}'"
    SUMMARY="${SUMMARY}  Target '${target}': ${extracted_pairs} read pairs\n"
done

echo ""
echo "=== Extraction Summary ==="
echo -e "$SUMMARY"
echo "Output files: ${OUTPUT_PREFIX}_*_R{1,2}.fq.gz"
echo "=== Done ==="
