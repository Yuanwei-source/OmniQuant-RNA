#!/usr/bin/env python3
"""Multi-species decontam ON vs OFF comparison.
Run from project root: python benchmark_results/scripts/decontam_comparison/compare_on_off.py
"""
import csv, math, os, sys
from collections import Counter

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(SCRIPT_DIR)))
os.chdir(PROJECT_ROOT)

SPECIES = {
    "drosophila_wolbachia": {
        "on": "benchmark_results/drosophila_wolbachia/runs/2026-05-27_drosophila_wolbachia_decontam-on/07.consensus_expression/Wolbachia_infected_vs_Wolbachia_free",
        "off": "benchmark_results/drosophila_wolbachia/runs/2026-05-27_drosophila_wolbachia_decontam-off/07.consensus_expression/Wolbachia_infected_vs_Wolbachia_free",
        "label": "Drosophila (Wolb_inf vs Free)"
    },
    "bombyx_mori": {
        "on": "benchmark_results/bombyx_mori/runs/2026-05-24_bombyx_mori_decontam-on/07.consensus_expression/Testis_vs_Ovary",
        "off": "benchmark_results/bombyx_mori/runs/2026-05-24_bombyx_mori_decontam-off/07.consensus_expression/Testis_vs_Ovary",
        "label": "Bombyx (Testis vs Ovary)"
    },
    "epicauta_diapause": {
        "on": "benchmark_results/epicauta_diapause/runs/2026-05-27_epicauta_diapause_decontam-on/07.consensus_expression/diapause_vs_non-diapause",
        "off": "benchmark_results/epicauta_diapause/runs/2026-05-27_epicauta_diapause_decontam-off/07.consensus_expression/diapause_vs_non-diapause",
        "label": "Epicauta (Diapause vs Non)"
    }
}

def load_consensus(path):
    genes = {}
    fpath = os.path.join(path, "consensus_results.tsv")
    with open(fpath) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            gid = row["gene_id_standard"]
            try: lfc = float(row["logFC_mean"]) if row["logFC_mean"] else None
            except: lfc = None
            genes[gid] = {"tier": row["tier"], "lfc": lfc}
    return genes

def analyze(label, on, off):
    tier_on = Counter(v["tier"] for v in on.values())
    tier_off = Counter(v["tier"] for v in off.values())
    common = set(on) & set(off)
    deg_on = {g for g in common if on[g]["tier"] in ("Tier_A", "Tier_B")}
    deg_off = {g for g in common if off[g]["tier"] in ("Tier_A", "Tier_B")}
    S, ONonly, OFFonly = deg_on & deg_off, deg_on - deg_off, deg_off - deg_on
    ta_both = sum(1 for g in common if on[g]["tier"] == "Tier_A" and off[g]["tier"] == "Tier_A")
    ta_on = tier_on.get("Tier_A", 0)
    ta_off = tier_off.get("Tier_A", 0)
    ratio = len(OFFonly) / len(ONonly) if len(ONonly) > 0 else float('inf')
    
    valid = [(on[g]["lfc"], off[g]["lfc"]) for g in common if on[g]["lfc"] and off[g]["lfc"]]
    if valid:
        xs, ys = [x for x,_ in valid], [y for _,y in valid]
        n = len(xs); mx, my = sum(xs)/n, sum(ys)/n
        num = sum((x-mx)*(y-my) for x,y in valid)
        den = math.sqrt(sum((x-mx)**2 for x in xs)) * math.sqrt(sum((y-my)**2 for y in ys))
        r = num/den if den else 0
    else:
        r, n = 0, 0

    delta_pct = (len(deg_off) - len(deg_on)) / len(deg_on) * 100 if len(deg_on) > 0 else 0
    jaccard = ta_both / (ta_on + ta_off - ta_both) if (ta_on + ta_off - ta_both) else 0

    print(f"\n{'='*65}")
    print(f" {label}")
    print(f"{'='*65}")
    print(f"  Universe={len(on)}")
    print(f"  ON  TA={ta_on} B={tier_on.get('Tier_B',0)} C={tier_on.get('Tier_C',0)} unc={tier_on.get('unclassified',0)}")
    print(f"  OFF TA={ta_off} B={tier_off.get('Tier_B',0)} C={tier_off.get('Tier_C',0)} unc={tier_off.get('unclassified',0)}")
    print(f"  DEG(A+B) ON={len(deg_on)} OFF={len(deg_off)} ({delta_pct:+.1f}%)")
    print(f"  Shared={len(S)} ON-only={len(ONonly)} OFF-only={len(OFFonly)} ratio={ratio:.1f}x")
    print(f"  TierA retention={100*ta_both/ta_on:.2f}% Jaccard={jaccard:.3f}")
    print(f"  logFC r={r:.4f} (n={n})")

    tiers = ["Tier_A", "Tier_B", "Tier_C", "unclassified"]
    cross = Counter()
    for g in common: cross[(on[g]["tier"], off[g]["tier"])] += 1
    print(f"\n  Tier transition (ON → OFF):")
    header = f"  {'ON\\OFF':<12}" + "".join(f"  {t:>8}" for t in tiers) + f"  {'sum':>6}"
    print(header)
    for t_on in tiers:
        rs = sum(cross[(t_on,to)] for to in tiers)
        row = f"  {t_on:<12}"
        for t_off in tiers:
            v = cross[(t_on,t_off)]
            row += f"  {v:>5}({100*v/rs:>3.0f}%)" if rs and v else f"  {v:>5}    "
        row += f"  {rs:>5}"
        print(row)

if __name__ == "__main__":
    for name, paths in SPECIES.items():
        on = load_consensus(paths["on"])
        off = load_consensus(paths["off"])
        analyze(paths["label"], on, off)
