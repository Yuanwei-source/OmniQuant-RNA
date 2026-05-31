#!/usr/bin/env python3
"""Run ablation analysis for any species. Takes species name as argument."""
import sys, os, csv, math
from collections import defaultdict

# Species configs
SPECIES = {
    "drosophila": {
        "base": "benchmark_results/drosophila_wolbachia/runs/2026-05-27_drosophila_wolbachia_decontam-on",
        "contrast": "Wolbachia_infected_vs_Wolbachia_free",
        "label": "Drosophila"
    },
    "bombyx": {
        "base": "benchmark_results/bombyx_mori/runs/2026-05-24_bombyx_mori_decontam-on",
        "contrast": "Testis_vs_Ovary",
        "label": "Bombyx"
    },
    "epicauta": {
        "base": "benchmark_results/epicauta_diapause/runs/2026-05-27_epicauta_diapause_decontam-on",
        "contrast": "diapause_vs_non-diapause",
        "label": "Epicauta"
    }
}

species_name = sys.argv[1] if len(sys.argv) > 1 else "drosophila"
cfg = SPECIES[species_name]
BASE = cfg["base"]
CONTRAST = cfg["contrast"]
LABEL = cfg["label"]

print(f"Running ablation for: {LABEL}")
print(f"  Base: {BASE}")
print(f"  Contrast: {CONTRAST}")

def load_consensus(path):
    rows = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            for k in ("support_n", "sign_consistency_n"):
                row[k] = int(row.get(k, 0) or 0)
            for k in ("best_rra_fdr", "best_cct_fdr", "logFC_CV", "logFC_mean"):
                v = row.get(k, None)
                try: row[k] = float(v) if v and v != "NA" and v != "NaN" else 99.0
                except (ValueError, TypeError): row[k] = 99.0
            rows.append(row)
    return rows

def load_dea(path):
    genes = {}
    if not os.path.exists(path):
        print(f"  WARNING: DEA file not found: {path}", file=sys.stderr)
        return genes
    with open(path) as f:
        for row in csv.DictReader(f):
            gid = row.get("gene_id_standard") or row.get("gene_id", "")
            if not gid: continue
            try: bm = float(row.get("baseMean", 0) or 0)
            except (ValueError, TypeError): bm = None
            try: lfc = abs(float(row.get("logFC", 0) or 0))
            except (ValueError, TypeError): lfc = 0
            try: padj = float(row.get("adj.P.Val", 1) or 1)
            except (ValueError, TypeError): padj = 1
            genes[gid] = {"baseMean": bm, "padj": padj, "abs_lfc": lfc}
    if not genes:
        print(f"  WARNING: No genes loaded from {path}", file=sys.stderr)
    return genes

def median(vals):
    if not vals: return 0
    s = sorted(vals); n = len(s)
    return s[n//2]

# Load
consensus_path = f"{BASE}/07.consensus_expression/{CONTRAST}/consensus_results.tsv"
dea_path = f"{BASE}/06.differential_expression/featurecounts/deseq2.{CONTRAST}.csv"

rows = load_consensus(consensus_path)
dea = load_dea(dea_path)
print(f"  Rows: {len(rows)}, DEA: {len(dea)}")

# Load per-quantifier sig status
quants = ["featurecounts", "stringtie", "salmon", "kallisto"]
sig_by_quant = defaultdict(lambda: defaultdict(bool))
for q in quants:
    import glob
    csvs = glob.glob(f"{BASE}/06.differential_expression/{q}/deseq2.*.csv")
    if csvs:
        for r in csv.DictReader(open(csvs[0])):
            gid = r.get("gene_id_standard") or r.get("gene_id", "")
            if gid:
                try:
                    sig = float(r.get("adj.P.Val", 1) or 1) < 0.05 and abs(float(r.get("logFC", 0) or 0)) > 1
                except (ValueError, TypeError):
                    sig = False
                sig_by_quant[gid][q] = sig

# Gene sets
any_signal = {g for g, sq in sig_by_quant.items() if sum(sq.values()) > 0}
full_ta = {r["gene_id_standard"] for r in rows if r["tier"] == "Tier_A"}
full_tb = {r["gene_id_standard"] for r in rows if r["tier"] == "Tier_B"}
vote_pass = {r["gene_id_standard"] for r in rows if r["support_n"] >= 3 and r["sign_consistency_n"] >= 3}
so_strict = {r["gene_id_standard"] for r in rows if r["support_n"] >= 4 and r["sign_consistency_n"] >= 4}
ta_no_rra = {r["gene_id_standard"] for r in rows if r["support_n"] >= 4 and r["sign_consistency_n"] >= 4
             and r["best_cct_fdr"] <= 0.05 and r["logFC_CV"] <= 1.0}

# Gates
fail_rra = {r["gene_id_standard"] for r in rows if r["tier"] != "Tier_A" and r["support_n"] >= 4 and r["sign_consistency_n"] >= 4
            and r["best_rra_fdr"] > 0.05 and r["best_cct_fdr"] <= 0.05 and r["logFC_CV"] <= 1.0
            and r.get("consensus_direction", "") not in ("mixed",)}
fail_cct = {r["gene_id_standard"] for r in rows if r["tier"] != "Tier_A" and r["support_n"] >= 4 and r["sign_consistency_n"] >= 4
            and r["best_cct_fdr"] > 0.05 and r["best_rra_fdr"] <= 0.05 and r["logFC_CV"] <= 1.0
            and r.get("consensus_direction", "") not in ("mixed",)}
fail_cv = {r["gene_id_standard"] for r in rows if r["tier"] != "Tier_A" and r["support_n"] >= 4 and r["sign_consistency_n"] >= 4
           and r["logFC_CV"] > 1.0 and r["best_rra_fdr"] <= 0.05 and r["best_cct_fdr"] <= 0.05
           and r.get("consensus_direction", "") not in ("mixed",)}
fail_rra_cct = {r["gene_id_standard"] for r in rows if r["tier"] != "Tier_A" and r["support_n"] >= 4 and r["sign_consistency_n"] >= 4
                and r["best_rra_fdr"] > 0.05 and r["best_cct_fdr"] > 0.05 and r["logFC_CV"] <= 1.0}

def pq(genes):
    exprs = [math.log10(dea.get(g, {}).get("baseMean", 0) + 1) for g in genes if dea.get(g, {}).get("baseMean")]
    lfcs = [dea.get(g, {}).get("abs_lfc", 0) for g in genes]
    fc_sig = sum(1 for g in genes if dea.get(g, {}).get("padj", 1) < 0.05)
    return {"n": len(genes), "med_expr": round(median(exprs), 2) if exprs else 0, "pct_fc": round(100*fc_sig/len(genes), 1) if genes else 0}

# OFF-only
on_tiers = defaultdict(set)
for r in rows: on_tiers[r["tier"]].add(r["gene_id_standard"])
deg_on = on_tiers["Tier_A"] | on_tiers["Tier_B"]
off_base = BASE.replace("-on", "-off")
off_rows = load_consensus(f"{off_base}/07.consensus_expression/{CONTRAST}/consensus_results.tsv")
off_tiers = defaultdict(set)
for r in off_rows: off_tiers[r["tier"]].add(r["gene_id_standard"])
deg_off = off_tiers["Tier_A"] | off_tiers["Tier_B"]
off_only = deg_off - deg_on

# Results
print(f"\n{'='*60}")
print(f" {LABEL}")
print(f"{'='*60}")
print(f"  Universe: {len(rows)}")
print(f"  Full Tier A: {len(full_ta)}  Tier B: {len(full_tb)}")
print(f"  Support-strict: {len(so_strict)}  Vote-pass: {len(vote_pass)}")
print(f"  Any signal: {len(any_signal)}")
print(f"\n  Gate failure (vote-pass candidate genome):")
print(f"    RRA only:    {len(fail_rra)}  {pq(fail_rra)}")
print(f"    CCT only:    {len(fail_cct)}  {pq(fail_cct)}")
print(f"    CV only:     {len(fail_cv)}  {pq(fail_cv)}")
print(f"    RRA+CCT:     {len(fail_rra_cct)}  {pq(fail_rra_cct)}")
print(f"\n  Support-strict vs Full TA:")
added = so_strict - full_ta
print(f"    Added by support-strict: {len(added)}  {pq(added)}")
removed = full_ta - so_strict
print(f"    Lost (Full TA not in support): {len(removed)}")
print(f"\n  OFF-only: {len(off_only)}  {pq(off_only)}")
print(f"\n  Quantifier overlap (Jaccard):")
qsets = {}
for q in quants:
    import glob
    csvs = glob.glob(f"{BASE}/06.differential_expression/{q}/deseq2.*.csv")
    if csvs:
        sigs = set()
        for r in csv.DictReader(open(csvs[0])):
            gid = r.get("gene_id_standard") or r.get("gene_id", "")
            if gid:
                try:
                    if float(r.get("adj.P.Val", 1) or 1) < 0.05 and abs(float(r.get("logFC", 0) or 0)) > 1:
                        sigs.add(gid)
                except (ValueError, TypeError):
                    pass
        qsets[q] = sigs
prefixes = {"featurecounts": "FC", "stringtie": "ST", "salmon": "SA", "kallisto": "KA"}
for q1 in quants:
    for q2 in quants:
        if prefixes[q1] < prefixes[q2]:
            j = len(qsets.get(q1, set()) & qsets.get(q2, set()))
            u = len(qsets.get(q1, set()) | qsets.get(q2, set()))
            print(f"    {prefixes[q1]} x {prefixes[q2]}: J={j/u:.3f}" if u else f"    {prefixes[q1]} x {prefixes[q2]}: N/A")
