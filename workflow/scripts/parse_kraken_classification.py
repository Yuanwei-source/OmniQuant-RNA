from __future__ import annotations

import argparse
import ast
import gzip
from collections import Counter, defaultdict
from pathlib import Path


def normalize_read_id(read_id: str) -> str:
    token = read_id.strip().split()[0]
    if token.startswith("@"):
        token = token[1:]
    if token.endswith("/1") or token.endswith("/2"):
        token = token[:-2]
    return token


def parse_taxid_list(value: str) -> set[str]:
    if not value or value in {"None", "null", "[]"}:
        return set()
    try:
        parsed = ast.literal_eval(value)
    except Exception:
        parsed = [part for part in value.replace(",", " ").split() if part]
    if isinstance(parsed, (list, tuple, set)):
        return {str(item) for item in parsed if str(item)}
    if parsed:
        return {str(parsed)}
    return set()


def load_parent_map(nodes_file: str) -> dict[str, str]:
    path = Path(nodes_file)
    if not path.exists():
        return {}

    parent_map: dict[str, str] = {}
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = [field.strip() for field in line.split("|")]
            if len(fields) < 2:
                continue
            taxid = fields[0]
            parent = fields[1]
            if taxid:
                parent_map[taxid] = parent
    return parent_map


def is_descendant(taxid: str, targets: set[str], parent_map: dict[str, str], cache: dict[tuple[str, tuple[str, ...]], bool]) -> bool:
    if not taxid or not targets:
        return False

    cache_key = (taxid, tuple(sorted(targets)))
    if cache_key in cache:
        return cache[cache_key]

    current = taxid
    seen = set()
    while current and current not in seen:
        if current in targets:
            cache[cache_key] = True
            return True
        seen.add(current)
        parent = parent_map.get(current)
        if not parent or parent == current:
            break
        current = parent

    cache[cache_key] = False
    return False


def write_id_gzip(output_path: str, ids: set[str]) -> None:
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(output_path, "wt") as handle:
        handle.write("# read_id\n")
        for read_id in sorted(ids):
            handle.write(f"{read_id}\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Parse Kraken2 assignments into nontarget and uncertain pair ID sets.")
    parser.add_argument("--kraken-output", required=True)
    parser.add_argument("--nodes-file", required=False)
    parser.add_argument("--taxonomy-nodes", required=False)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--nontarget-taxids", default="[]")
    parser.add_argument("--uncertain-taxids", default="[]")
    parser.add_argument("--ecological-as-uncertain", default="True")
    parser.add_argument("--output-nontarget", required=True)
    parser.add_argument("--output-uncertain", required=True)
    parser.add_argument("--output-stats", required=True)
    args = parser.parse_args()

    nontarget_taxids = parse_taxid_list(args.nontarget_taxids)
    uncertain_taxids = parse_taxid_list(args.uncertain_taxids)
    parent_map = load_parent_map(args.taxonomy_nodes or "")
    descendant_cache: dict[tuple[str, tuple[str, ...]], bool] = {}

    pair_flags: dict[str, dict[str, bool]] = defaultdict(lambda: {
        "nontarget": False,
        "uncertain": False,
        "classified": False,
        "unclassified": False,
    })
    counters = Counter()

    kraken_path = Path(args.kraken_output)
    if kraken_path.exists() and kraken_path.stat().st_size > 0:
        with kraken_path.open("r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue
                status = fields[0]
                read_id = normalize_read_id(fields[1])
                taxid = fields[2].strip()
                if not read_id:
                    continue

                if status == "C":
                    pair_flags[read_id]["classified"] = True
                    counters["classified_lines"] += 1
                    if is_descendant(taxid, nontarget_taxids, parent_map, descendant_cache):
                        pair_flags[read_id]["nontarget"] = True
                    elif is_descendant(taxid, uncertain_taxids, parent_map, descendant_cache):
                        pair_flags[read_id]["uncertain"] = True
                else:
                    pair_flags[read_id]["unclassified"] = True
                    counters["unclassified_lines"] += 1

    nontarget_ids: set[str] = set()
    uncertain_ids: set[str] = set()
    ecological_as_uncertain = str(args.ecological_as_uncertain).lower() == "true"

    for read_id, flags in pair_flags.items():
        if flags["nontarget"]:
            nontarget_ids.add(read_id)
        elif flags["uncertain"] and ecological_as_uncertain:
            uncertain_ids.add(read_id)

    write_id_gzip(args.output_nontarget, nontarget_ids)
    write_id_gzip(args.output_uncertain, uncertain_ids)

    with open(args.output_stats, "w", encoding="utf-8") as handle:
        handle.write("sample\tstage\tclassified_pairs\tunclassified_pairs\tnontarget_pairs\tuncertain_pairs\n")
        handle.write(
            "{sample}\tclassification\t{classified}\t{unclassified}\t{nontarget}\t{uncertain}\n".format(
                sample=args.sample,
                classified=sum(1 for flags in pair_flags.values() if flags["classified"]),
                unclassified=sum(1 for flags in pair_flags.values() if flags["unclassified"]),
                nontarget=len(nontarget_ids),
                uncertain=len(uncertain_ids),
            )
        )


if __name__ == "__main__":
    main()