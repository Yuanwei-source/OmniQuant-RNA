from collections import Counter
from pathlib import Path
import logging

from xopen import xopen


def setup_logging(log_path: str) -> None:
    Path(log_path).parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        filename=log_path,
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )


def normalize_read_id(header: str) -> str:
    token = header.strip().split()[0]
    if token.startswith("@"):
        token = token[1:]
    if token.endswith("/1") or token.endswith("/2"):
        token = token[:-2]
    return token


def load_id_set(path: str) -> set[str]:
    read_ids: set[str] = set()
    with xopen(path, "rt") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            read_id = normalize_read_id(line.split("\t", 1)[0])
            if read_id:
                read_ids.add(read_id)
    return read_ids


def fastq_pair_iter(r1_path: str, r2_path: str):
    with xopen(r1_path, "rt") as r1_handle, xopen(r2_path, "rt") as r2_handle:
        while True:
            record_r1 = [r1_handle.readline() for _ in range(4)]
            record_r2 = [r2_handle.readline() for _ in range(4)]

            if not any(record_r1) and not any(record_r2):
                break

            if any(line == "" for line in record_r1) or any(line == "" for line in record_r2):
                raise ValueError("FASTQ files ended unevenly; paired-end integrity is broken before decontamination.")

            read_id_r1 = normalize_read_id(record_r1[0])
            read_id_r2 = normalize_read_id(record_r2[0])
            if read_id_r1 != read_id_r2:
                raise ValueError(
                    f"FASTQ pair ID mismatch detected: '{read_id_r1}' vs '{read_id_r2}'. "
                    "The input paired-end files are not synchronized."
                )

            yield read_id_r1, record_r1, record_r2


def decide_fate(
    read_id: str,
    mode: str,
    retain_unclassified: bool,
    tech_ids: set[str],
    host_ids: set[str],
    ercc_ids: set[str],
    nontarget_ids: set[str],
    uncertain_ids: set[str],
) -> str:
    if read_id in ercc_ids:
        return "Retained_ERCC"
    if read_id in host_ids:
        return "Retained_Host"

    if mode == "report_only":
        return "Retained_Uncertain"

    if read_id in tech_ids:
        return "Removed_Tech"
    if read_id in nontarget_ids:
        return "Removed_NonTarget"
    if read_id in uncertain_ids or retain_unclassified:
        return "Retained_Uncertain"

    return "Removed_NonTarget"


def write_record_pair(handle_r1, handle_r2, record_r1, record_r2) -> None:
    handle_r1.writelines(record_r1)
    handle_r2.writelines(record_r2)


def main() -> None:
    log_path = str(snakemake.log[0])
    setup_logging(log_path)

    mode = str(snakemake.params.mode)
    retain_unclassified = bool(snakemake.params.retain_unclassified)
    keep_audit_fastq = bool(snakemake.params.keep_audit_fastq)

    tech_ids = load_id_set(str(snakemake.input.tech_ids))
    host_ids = load_id_set(str(snakemake.input.host_ids))
    ercc_ids = load_id_set(str(snakemake.input.ercc_ids))
    nontarget_ids = load_id_set(str(snakemake.input.nontarget_ids))
    uncertain_ids = load_id_set(str(snakemake.input.uncertain_ids))

    logging.info(
        "Loaded evidence ID sets for sample %s: tech=%d host=%d ercc=%d nontarget=%d uncertain=%d",
        snakemake.wildcards.sample,
        len(tech_ids),
        len(host_ids),
        len(ercc_ids),
        len(nontarget_ids),
        len(uncertain_ids),
    )

    output_paths = [
        str(snakemake.output.clean_r1),
        str(snakemake.output.clean_r2),
        str(snakemake.output.uncertain_r1),
        str(snakemake.output.uncertain_r2),
        str(snakemake.output.removed_r1),
        str(snakemake.output.removed_r2),
        str(snakemake.output.summary),
    ]
    for output_path in output_paths:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    fate_counts: Counter[str] = Counter()
    total_pairs = 0

    with xopen(str(snakemake.output.clean_r1), "wt") as clean_r1_handle, \
        xopen(str(snakemake.output.clean_r2), "wt") as clean_r2_handle, \
        xopen(str(snakemake.output.uncertain_r1), "wt") as uncertain_r1_handle, \
        xopen(str(snakemake.output.uncertain_r2), "wt") as uncertain_r2_handle, \
        xopen(str(snakemake.output.removed_r1), "wt") as removed_r1_handle, \
        xopen(str(snakemake.output.removed_r2), "wt") as removed_r2_handle:

        for read_id, record_r1, record_r2 in fastq_pair_iter(str(snakemake.input.r1), str(snakemake.input.r2)):
            total_pairs += 1
            fate = decide_fate(
                read_id=read_id,
                mode=mode,
                retain_unclassified=retain_unclassified,
                tech_ids=tech_ids,
                host_ids=host_ids,
                ercc_ids=ercc_ids,
                nontarget_ids=nontarget_ids,
                uncertain_ids=uncertain_ids,
            )
            fate_counts[fate] += 1

            if fate in {"Retained_Host", "Retained_ERCC", "Retained_Uncertain"}:
                write_record_pair(clean_r1_handle, clean_r2_handle, record_r1, record_r2)
                if keep_audit_fastq and fate == "Retained_Uncertain":
                    write_record_pair(uncertain_r1_handle, uncertain_r2_handle, record_r1, record_r2)
            else:
                if keep_audit_fastq:
                    write_record_pair(removed_r1_handle, removed_r2_handle, record_r1, record_r2)

    with open(str(snakemake.output.summary), "w", encoding="utf-8") as handle:
        handle.write(
            "Sample\tRetained_Host\tRetained_ERCC\tRetained_Uncertain\tRemoved_Tech\tRemoved_NonTarget\n"
        )
        handle.write(
            "{sample}\t{host}\t{ercc}\t{uncertain}\t{tech}\t{nontarget}\n".format(
                sample=snakemake.wildcards.sample,
                host=fate_counts.get("Retained_Host", 0),
                ercc=fate_counts.get("Retained_ERCC", 0),
                uncertain=fate_counts.get("Retained_Uncertain", 0),
                tech=fate_counts.get("Removed_Tech", 0),
                nontarget=fate_counts.get("Removed_NonTarget", 0),
            )
        )

    logging.info(
        "Finished pair-aware decision for sample %s: total_pairs=%d host=%d ercc=%d uncertain=%d tech_removed=%d nontarget_removed=%d",
        snakemake.wildcards.sample,
        total_pairs,
        fate_counts.get("Retained_Host", 0),
        fate_counts.get("Retained_ERCC", 0),
        fate_counts.get("Retained_Uncertain", 0),
        fate_counts.get("Removed_Tech", 0),
        fate_counts.get("Removed_NonTarget", 0),
    )


main()