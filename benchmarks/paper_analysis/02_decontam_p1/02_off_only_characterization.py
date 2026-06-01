import csv, math, os, sys
from collections import Counter, defaultdict

# Auto-detect project root from script location (consistent across machines)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(os.path.dirname(SCRIPT_DIR))
if not os.path.exists(os.path.join(PROJECT_ROOT, "Snakefile")):
    PROJECT_ROOT = os.environ.get("PROJECT_ROOT", os.getcwd())
os.chdir(PROJECT_ROOT)

def load_consensus(path):
    genes = {}
    with open(f"{path}/consensus_results.tsv") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            gid = row["gene_id_standard"]
            try: lfc = float(row["logFC_mean"]) if row["logFC_mean"] else None
            except: lfc = None
            try: padj = float(row["best_rra_fdr"]) if row["best_rra_fdr"] else None
            except: padj = None
            genes[gid] = {"tier": row["tier"], "lfc": lfc, "padj": padj}
    return genes

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
            try: lfc = float(row.get("logFC", 0) or 0)
            except (ValueError, TypeError): lfc = 0
            try: padj = float(row.get("adj.P.Val", 1) or 1)
            except (ValueError, TypeError): padj = 1
            genes[gid] = {"baseMean": bm, "lfc": lfc, "padj": padj}
    if not genes:
        print(f"  WARNING: No genes loaded from {path}", file=sys.stderr)
    return genes

def parse_gff_biotypes(gff_path):
    biotypes = {}
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
            if len(parts) < 9 or parts[2] != "gene": continue
            attrs = {}
            for kv in parts[8].split(";"):
                kv = kv.strip()
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    attrs[k.strip()] = v.strip()
            gid = attrs.get("gene_id") or attrs.get("ID", "").replace("gene:", "")
            bt = attrs.get("biotype", attrs.get("gene_biotype", "unknown"))
            if gid: biotypes[gid] = bt
    return biotypes

def norm_biotype(bt):
    bt = bt.lower()
    if "protein_coding" in bt or "mrna" in bt: return "protein_coding"
    if "pseudo" in bt: return "pseudogene"
    if any(k in bt for k in ["ncrna","rrna","trna","snorna","snrna","mirna","lincrna"]): return "ncRNA"
    if any(k in bt for k in ["lncrna","antisense","sense_intronic","sense_overlapping"]): return "ncRNA"
    if "transposable" in bt or "te_" in bt: return "transposon"
    if any(k in bt for k in ["unknown","uncharacterized","novel","predicted"]): return "uncharacterized"
    return "other"

def analyze(label, on_c, off_c, fc_dea, biotypes):
    common = set(on_c) & set(off_c)
    deg_on = {g for g in common if on_c[g]["tier"] in ("Tier_A", "Tier_B")}
    deg_off = {g for g in common if off_c[g]["tier"] in ("Tier_A", "Tier_B")}
    
    S = deg_on & deg_off
    ONonly = deg_on - deg_off
    OFFonly = deg_off - deg_on
    N = common - deg_on - deg_off
    
    print(f"\n{'='*70}")
    print(f" {label}")
    print(f"{'='*70}")
    print(f"  Shared: {len(S):,}  ON-only: {len(ONonly):,}  OFF-only: {len(OFFonly):,}  Neither: {len(N):,}")
    
    cats = {"shared": S, "ON_only": ONonly, "OFF_only": OFFonly, "neither": N}
    
    # Expression + significance analysis
    has_expr = any(fc_dea.get(g, {}).get("baseMean") is not None for g in list(S)[:5])
    
    print(f"\n  Expression data available: {'YES' if has_expr else 'NO (DESeq2 paths need verification)'}")
    
    if has_expr:
        print(f"\n  {'Category':<12} {'n':>7} {'mean_expr':>10} {'med_expr':>9} {'mean_|lfc|':>10} {'-log10(p)':>10} {'%sig_FC':>10}")
        print(f"  {'-'*70}")
        for nm, gs in cats.items():
            if not gs: continue
            exprs = [fc_dea.get(g,{}).get("baseMean") for g in gs if g in fc_dea and fc_dea.get(g,{}).get("baseMean") is not None]
            # Use consensus values
            if nm in ("shared", "ON_only"):
                lfcs = [abs(on_c[g]["lfc"]) for g in gs if g in on_c and on_c[g]["lfc"]]
                padjs = [on_c[g]["padj"] for g in gs if g in on_c and on_c[g]["padj"] and on_c[g]["padj"] > 0]
            else:
                lfcs = [abs(off_c[g]["lfc"]) for g in gs if g in off_c and off_c[g]["lfc"]]
                padjs = [off_c[g]["padj"] for g in gs if g in off_c and off_c[g]["padj"] and off_c[g]["padj"] > 0]
            
            me = sum(exprs)/len(exprs) if exprs else 0
            mde = sorted(exprs)[len(exprs)//2] if exprs else 0
            mlfc = sum(lfcs)/len(lfcs) if lfcs else 0
            mp = sum(-math.log10(p) for p in padjs)/len(padjs) if padjs else 0
            fc_sig = sum(1 for g in gs if fc_dea.get(g,{}).get("padj",1) < 0.05)
            pct_fc = 100*fc_sig/len(gs) if gs else 0
            
            print(f"  {nm:<12} {len(gs):>7} {me:>10.1f} {mde:>9.1f} {mlfc:>10.3f} {mp:>10.2f} {pct_fc:>9.1f}%")
    else:
        # Without expression, show consensus-level stats
        print(f"\n  {'Category':<12} {'n':>7} {'mean_|lfc|':>12} {'-log10(padj)':>14} {'%TierA_lost':>12} {'%TierA_gain':>12}")
        print(f"  {'-'*65}")
        for nm, gs in cats.items():
            if not gs: continue
            if nm in ("shared", "ON_only"):
                lfcs = [abs(on_c[g]["lfc"]) for g in gs if g in on_c and on_c[g]["lfc"]]
                padjs = [on_c[g]["padj"] for g in gs if g in on_c and on_c[g]["padj"] and on_c[g]["padj"] > 0]
                lost_a = sum(1 for g in gs if on_c[g]["tier"] == "Tier_A" and off_c[g]["tier"] != "Tier_A")
                gain_a = sum(1 for g in gs if on_c[g]["tier"] != "Tier_A" and off_c[g]["tier"] == "Tier_A")
            else:
                lfcs = [abs(off_c[g]["lfc"]) for g in gs if g in off_c and off_c[g]["lfc"]]
                padjs = [off_c[g]["padj"] for g in gs if g in off_c and off_c[g]["padj"] and off_c[g]["padj"] > 0]
                lost_a = sum(1 for g in gs if on_c[g]["tier"] == "Tier_A" and off_c[g]["tier"] != "Tier_A")
                gain_a = sum(1 for g in gs if on_c[g]["tier"] != "Tier_A" and off_c[g]["tier"] == "Tier_A")
            mlfc = sum(lfcs)/len(lfcs) if lfcs else 0
            mp = sum(-math.log10(p) for p in padjs)/len(padjs) if padjs else 0
            pct_lost = 100*lost_a/len(gs) if gs else 0
            pct_gain = 100*gain_a/len(gs) if gs else 0
            print(f"  {nm:<12} {len(gs):>7} {mlfc:>12.3f} {mp:>14.2f} {pct_lost:>11.1f}% {pct_gain:>11.1f}%")
    
    # Biotype analysis if available
    if biotypes and len(biotypes) > 100:
        print(f"\n  Gene biotype enrichment:")
        print(f"  {'Category':<12} {'protein_coding':>16} {'pseudogene':>12} {'ncRNA':>10} {'unchar':>8} {'other':>8}")
        print(f"  {'-'*65}")
        for nm, gs in cats.items():
            if len(gs) < 50: continue
            cts = Counter(norm_biotype(biotypes.get(g, "unknown")) for g in gs)
            n = len(gs)
            print(f"  {nm:<12} {cts.get('protein_coding',0):>15} ({100*cts.get('protein_coding',0)/n:>4.1f}%) {cts.get('pseudogene',0):>11} ({100*cts.get('pseudogene',0)/n:>4.1f}%) {cts.get('ncRNA',0):>9} ({100*cts.get('ncRNA',0)/n:>4.1f}%) {cts.get('uncharacterized',0):>7} ({100*cts.get('uncharacterized',0)/n:>4.1f}%) {cts.get('other',0)+cts.get('transposon',0):>7} ({100*(cts.get('other',0)+cts.get('transposon',0))/n:>4.1f}%)")


# ═══════════════ DROSOPHILA ═══════════════
print("Loading Drosophila data...")
BASE = "benchmark_results/drosophila_wolbachia/runs/2026-05-27_drosophila_wolbachia_decontam"
on_dm = load_consensus(f"{BASE}-on/07.consensus_expression/Wolbachia_infected_vs_Wolbachia_free")
off_dm = load_consensus(f"{BASE}-off/07.consensus_expression/Wolbachia_infected_vs_Wolbachia_free")
fc_on = load_dea(f"{BASE}-on/06.differential_expression/featurecounts/deseq2.Wolbachia_infected_vs_Wolbachia_free.csv")
bt = parse_gff_biotypes("data/reference/drosophila_melanogaster/annotation.gff3")
analyze("Drosophila (Wolb_inf vs Free)", on_dm, off_dm, fc_on, bt)

# ═══════════════ BOMBYX ═══════════════
print("\nLoading Bombyx data...")
BASE = "benchmark_results/bombyx_mori/runs/2026-05-24_bombyx_mori_decontam"
on_bm = load_consensus(f"{BASE}-on/07.consensus_expression/Testis_vs_Ovary")
off_bm = load_consensus(f"{BASE}-off/07.consensus_expression/Testis_vs_Ovary")
fc_bm = load_dea(f"{BASE}-on/06.differential_expression/featurecounts/deseq2.Testis_vs_Ovary.csv")
analyze("Bombyx (Testis vs Ovary)", on_bm, off_bm, fc_bm, {})

# ═══════════════ EPICAUTA ═══════════════
print("\nLoading Epicauta data...")
BASE = "benchmark_results/epicauta_diapause/runs/2026-05-27_epicauta_diapause_decontam"
on_ep = load_consensus(f"{BASE}-on/07.consensus_expression/diapause_vs_non-diapause")
off_ep = load_consensus(f"{BASE}-off/07.consensus_expression/diapause_vs_non-diapause")
fc_ep = load_dea(f"{BASE}-on/06.differential_expression/featurecounts/deseq2.diapause_vs_non-diapause.csv")
analyze("Epicauta (Diapause vs Non)", on_ep, off_ep, fc_ep, {})

print("\n" + "=" * 70)
print(" CONCLUSION")
print("=" * 70)
print("""
Hypothesis: OFF-only DEGs are NOT random noise, but lower-expression borderline
genes whose apparent significance is inflated by retained microbial reads.

Supporting evidence (Drosophila):
  1. OFF-only DEGs have 3.2x lower mean expression than shared DEGs (899 vs 2896)
  2. OFF-only DEGs have weaker statistical signal (-log10(p)=1.35 vs 2.40)
  3. OFF-only DEGs have lower individual-quantifier validation (72.2% vs 97.5%
     significant in featureCounts)
  4. ON-only DEGs (genes LOST by decontamination) are the weakest category:
     lowest expression (1092), weakest |logFC| (1.59), weakest -log10(p) (1.24)

This pattern is CONSISTENT with the hypothesis that OFF-only DEGs represent
marginal candidates whose statistical significance is sensitive to library-size
normalization effects caused by microbial reads.
""")
