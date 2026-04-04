import csv
from pathlib import Path


DEFAULT_PRIORITY_TARGETS = {
    "wolbachia": "953",
    "virus": "10239",
    "fungi": "4751",
    "environmental_bacteria": "2",
}

COMPOSITION_CATEGORIES = [
    "Technical",
    "Bacteria",
    "Fungi",
    "Viruses",
    "Other_Classified",
    "Unclassified",
]


def read_sample_sheet(path):
    with open(path, "r", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sample_from_suffix(path_str, suffix):
    name = Path(path_str).name
    if not name.endswith(suffix):
        raise ValueError(f"Unexpected file name for suffix {suffix}: {name}")
    return name[: -len(suffix)]


def parse_decision_summary(path):
    expected_header = [
        "sample",
        "Rescued_Host",
        "Rescued_ERCC",
        "Flagged_Uncertain",
        "Removed_Tech",
        "Removed_NonTarget",
    ]

    with open(path, "r", encoding="utf-8") as handle:
        rows = [line.strip().split("\t") for line in handle if line.strip()]

    if not rows:
        raise ValueError(f"Empty decision summary: {path}")

    header_like = rows[0][0].lower() == "sample"
    data_row = rows[1] if header_like and len(rows) > 1 else rows[0]
    if len(data_row) != len(expected_header):
        raise ValueError(f"Unexpected decision summary shape in {path}: {data_row}")

    record = dict(zip(expected_header, data_row))
    return {
        "sample": record["sample"],
        "rescued_host_pairs": int(record["Rescued_Host"]),
        "rescued_ercc_pairs": int(record["Rescued_ERCC"]),
        "flagged_uncertain_pairs": int(record["Flagged_Uncertain"]),
        "removed_tech_pairs": int(record["Removed_Tech"]),
        "removed_nontarget_pairs": int(record["Removed_NonTarget"]),
    }


def parse_classification_stats(path):
    with open(path, "r", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    if not rows:
        return {
            "classified_pairs": 0,
            "unclassified_pairs": 0,
            "nontarget_pairs": 0,
            "uncertain_pairs": 0,
        }

    row = rows[0]
    return {
        "classified_pairs": int(row.get("classified_pairs", 0) or 0),
        "unclassified_pairs": int(row.get("unclassified_pairs", 0) or 0),
        "nontarget_pairs": int(row.get("nontarget_pairs", 0) or 0),
        "uncertain_pairs": int(row.get("uncertain_pairs", 0) or 0),
    }


def parse_kraken_report(path):
    rows = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 6:
                continue
            rows.append(
                {
                    "percent": float(parts[0]),
                    "clade_reads": int(float(parts[1])),
                    "direct_reads": int(float(parts[2])),
                    "rank_code": parts[-3].strip(),
                    "taxid": parts[-2].strip(),
                    "name": parts[-1].strip(),
                }
            )
    return rows


def safe_fraction(numerator, denominator):
    if not denominator:
        return 0.0
    return numerator / denominator


def classify_strength(fraction):
    if fraction >= 0.05:
        return "high"
    if fraction >= 0.01:
        return "medium"
    if fraction > 0:
        return "low"
    return "absent"


def top_taxa_summary(report_rows, top_n):
    filtered = [
        row for row in report_rows
        if row["taxid"] not in {"0", "1"}
        and row["name"].lower() not in {"root", "unclassified"}
        and row["clade_reads"] > 0
    ]
    filtered.sort(key=lambda row: (-row["clade_reads"], row["name"]))
    return "; ".join(f"{row['name']}({row['clade_reads']})" for row in filtered[:top_n])


sample_rows = read_sample_sheet(snakemake.input.sample_file)
sample_order = [row["sample"] for row in sample_rows]
group_map = {row["sample"]: row.get("group", "") for row in sample_rows}

decision_map = {
    sample_from_suffix(path, "_decision_summary.tsv"): parse_decision_summary(path)
    for path in snakemake.input.decision_summaries
}
classification_map = {
    sample_from_suffix(path, "_classification_stats.tsv"): parse_classification_stats(path)
    for path in snakemake.input.classification_stats
}
report_map = {
    sample_from_suffix(path, ".kraken.report.tsv"): parse_kraken_report(path)
    for path in snakemake.input.reports
}

clues_config = snakemake.config.get("decontam", {}).get("clues", {})
target_taxids = {
    key: str(value)
    for key, value in clues_config.get("priority_targets", DEFAULT_PRIORITY_TARGETS).items()
}
top_n = int(clues_config.get("top_taxa_n", 3))

burden_fields = [
    "sample",
    "group",
    "total_pairs",
    "rescued_host_pairs",
    "rescued_ercc_pairs",
    "non_host_pairs",
    "non_host_fraction",
    "technical_fraction",
    "nontarget_fraction",
    "flagged_uncertain_fraction",
    "classified_pairs",
    "unclassified_pairs",
    "bacteria_pairs",
    "bacteria_fraction",
    "fungi_pairs",
    "fungi_fraction",
    "virus_pairs",
    "virus_fraction",
    "top_taxa_summary",
]

target_fields = ["sample", "group"]
for label in ["wolbachia", "virus", "fungi", "environmental_bacteria"]:
    target_fields.extend(
        [
            f"{label}_pairs",
            f"{label}_fraction",
            f"{label}_presence",
            f"{label}_strength",
        ]
    )
target_fields.extend(["dominant_nonhost_taxon", "dominant_nonhost_pairs"])

burden_rows = []
target_rows = []
composition_rows = []

for sample in sample_order:
    decision = decision_map[sample]
    classification = classification_map.get(sample, {})
    report_rows = report_map.get(sample, [])
    report_by_taxid = {row["taxid"]: row for row in report_rows}

    total_pairs = (
        decision["rescued_host_pairs"]
        + decision["rescued_ercc_pairs"]
        + decision["flagged_uncertain_pairs"]
        + decision["removed_tech_pairs"]
        + decision["removed_nontarget_pairs"]
    )
    non_host_pairs = total_pairs - decision["rescued_host_pairs"] - decision["rescued_ercc_pairs"]

    bacteria_pairs = report_by_taxid.get("2", {}).get("clade_reads", 0)
    fungi_pairs = report_by_taxid.get("4751", {}).get("clade_reads", 0)
    virus_pairs = report_by_taxid.get("10239", {}).get("clade_reads", 0)
    classified_pairs = classification.get("classified_pairs", 0)
    unclassified_pairs = classification.get("unclassified_pairs", 0)
    other_classified_pairs = max(classified_pairs - bacteria_pairs - fungi_pairs - virus_pairs, 0)

    burden_rows.append(
        {
            "sample": sample,
            "group": group_map.get(sample, ""),
            "total_pairs": total_pairs,
            "rescued_host_pairs": decision["rescued_host_pairs"],
            "rescued_ercc_pairs": decision["rescued_ercc_pairs"],
            "non_host_pairs": non_host_pairs,
            "non_host_fraction": f"{safe_fraction(non_host_pairs, total_pairs):.6f}",
            "technical_fraction": f"{safe_fraction(decision['removed_tech_pairs'], total_pairs):.6f}",
            "nontarget_fraction": f"{safe_fraction(decision['removed_nontarget_pairs'], total_pairs):.6f}",
            "flagged_uncertain_fraction": f"{safe_fraction(decision['flagged_uncertain_pairs'], total_pairs):.6f}",
            "classified_pairs": classified_pairs,
            "unclassified_pairs": unclassified_pairs,
            "bacteria_pairs": bacteria_pairs,
            "bacteria_fraction": f"{safe_fraction(bacteria_pairs, total_pairs):.6f}",
            "fungi_pairs": fungi_pairs,
            "fungi_fraction": f"{safe_fraction(fungi_pairs, total_pairs):.6f}",
            "virus_pairs": virus_pairs,
            "virus_fraction": f"{safe_fraction(virus_pairs, total_pairs):.6f}",
            "top_taxa_summary": top_taxa_summary(report_rows, top_n),
        }
    )

    target_row = {
        "sample": sample,
        "group": group_map.get(sample, ""),
    }
    for label in ["wolbachia", "virus", "fungi", "environmental_bacteria"]:
        target_taxid = target_taxids[label]
        pairs = report_by_taxid.get(target_taxid, {}).get("clade_reads", 0)
        fraction = safe_fraction(pairs, total_pairs)
        target_row[f"{label}_pairs"] = pairs
        target_row[f"{label}_fraction"] = f"{fraction:.6f}"
        target_row[f"{label}_presence"] = "yes" if pairs > 0 else "no"
        target_row[f"{label}_strength"] = classify_strength(fraction)

    top_taxa = [
        row for row in report_rows
        if row["taxid"] not in {"0", "1"} and row["name"].lower() not in {"root", "unclassified"}
    ]
    top_taxa.sort(key=lambda row: (-row["clade_reads"], row["name"]))
    top_hit = top_taxa[0] if top_taxa else {"name": "none", "clade_reads": 0}
    target_row["dominant_nonhost_taxon"] = top_hit["name"]
    target_row["dominant_nonhost_pairs"] = top_hit["clade_reads"]
    target_rows.append(target_row)

    composition_values = {
        "Technical": decision["removed_tech_pairs"],
        "Bacteria": bacteria_pairs,
        "Fungi": fungi_pairs,
        "Viruses": virus_pairs,
        "Other_Classified": other_classified_pairs,
        "Unclassified": unclassified_pairs,
    }
    for category in COMPOSITION_CATEGORIES:
        pairs = composition_values[category]
        composition_rows.append(
            {
                "sample": sample,
                "group": group_map.get(sample, ""),
                "category": category,
                "pairs": pairs,
                "fraction_total": f"{safe_fraction(pairs, total_pairs):.6f}",
                "fraction_non_host": f"{safe_fraction(pairs, non_host_pairs):.6f}",
            }
        )

for output_path in snakemake.output:
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

with open(snakemake.output.burden, "w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=burden_fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(burden_rows)

with open(snakemake.output.targets, "w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=target_fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(target_rows)

with open(snakemake.output.composition, "w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=["sample", "group", "category", "pairs", "fraction_total", "fraction_non_host"],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(composition_rows)