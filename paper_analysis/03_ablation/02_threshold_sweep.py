import csv
from collections import Counter
import sys

RUN = "benchmark_results/drosophila_wolbachia/runs/2026-05-27_drosophila_wolbachia_decontam-on/07.consensus_expression/Wolbachia_infected_vs_Wolbachia_free"

# Read tier_diagnostics
rows = []
with open(f"{RUN}/tier_diagnostics.tsv") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        rows.append(row)

total = len(rows)
tier_a = sum(1 for r in rows if r["tier"] == "Tier_A")
tier_b = sum(1 for r in rows if r["tier"] == "Tier_B")
tier_c = sum(1 for r in rows if r["tier"] == "Tier_C")
unclass = sum(1 for r in rows if r["tier"] == "unclassified")

print(f"=== 第1部分: Tier A 失败原因分析 ===")
print(f"总基因数: {total}")
print(f"Tier A: {tier_a} ({100*tier_a/total:.1f}%)")
print(f"未达 Tier A: {total-tier_a} ({100*(total-tier_a)/total:.1f}%)\n")

# Parse tier_blocker_a
blocker_a = Counter()
for r in rows:
    if r["tier"] != "Tier_A":
        blocker = r["tier_blocker_a"]
        blocker_a[blocker] += 1

print("Tier_A 阻断因素:")
for k, v in blocker_a.most_common():
    pct = 100*v/(total-tier_a)
    print(f"  {k:<60} {v:>6} ({pct:5.1f}%)")

# Group by category
no_signal = blocker_a.get("no_signal", 0)
rra_related = sum(v for k,v in blocker_a.items() if "rra_fdr" in k)
cct_related = sum(v for k,v in blocker_a.items() if "cct_fdr" in k)
support_related = sum(v for k,v in blocker_a.items() if "support" in k)
mixed = blocker_a.get("mixed_direction", 0)
cv_related = sum(v for k,v in blocker_a.items() if "logFC_CV" in k)

print(f"\n分类汇总:")
print(f"  无任何定量器信号: {no_signal:>6} ({100*no_signal/(total-tier_a):.1f}%)")
print(f"  RRA 相关失败:    {rra_related:>6} ({100*rra_related/(total-tier_a):.1f}%)")
print(f"  CCT 相关失败:    {cct_related:>6} ({100*cct_related/(total-tier_a):.1f}%)")
print(f"  支持度/一致性:    {support_related:>6} ({100*support_related/(total-tier_a):.1f}%)")
print(f"  方向冲突:        {mixed:>6} ({100*mixed/(total-tier_a):.1f}%)")
print(f"  logFC_CV 失败:   {cv_related:>6} ({100*cv_related/(total-tier_a):.1f}%)")

# Purely statistical (only RRA/CCT, no support/consistency)
pure_stat = sum(v for k,v in blocker_a.items() if k in ("rra_fdr", "cct_fdr", "rra_fdr;cct_fdr"))
print(f"\n纯统计门控失败（已通过支持度+一致性，仅卡在 P 值）: {pure_stat} ({100*pure_stat/(total-tier_a):.1f}%)")

# Tier A that pass all statistical gates but fail logFC_CV
cv_only = sum(v for k,v in blocker_a.items() if "logFC_CV" in k and "rra_fdr" not in k and "cct_fdr" not in k and "support" not in k)
print(f"独立因 logFC_CV 失败: {cv_only}")

print(f"\n=== 第2部分: Tier B 门控敏感性 ===")
# Genes currently in Tier B
tier_b_genes = [r for r in rows if r["tier"] == "Tier_B"]
print(f"当前 Tier B 基因数: {len(tier_b_genes)}")

# For tier B genes, count what would happen if each gate were removed
# We need: support_n, sign_consistency_n, best_rra_fdr, best_cct_fdr, logFC_CV
# Tier B thresholds: support>=3, sign_consistency>=3, rra_fdr<0.10, cct_fdr<0.10, logFC_CV<1.25

# Scenario: without RRA (only CCT applies)
without_rra = 0
skipped_rra = 0
for r in tier_b_genes:
    try:
        support = int(r["support_n"])
        sign = int(r["sign_consistency_n"])
        cct = float(r["best_cct_fdr"])
        cv = float(r["logFC_CV"])
        if support >= 3 and sign >= 3 and cct < 0.10 and cv < 1.25:
            without_rra += 1
    except (ValueError, TypeError, KeyError):
        skipped_rra += 1

# Scenario: without CCT (only RRA applies)
without_cct = 0
skipped_cct = 0
for r in tier_b_genes:
    try:
        support = int(r["support_n"])
        sign = int(r["sign_consistency_n"])
        rra = float(r["best_rra_fdr"])
        cv = float(r["logFC_CV"])
        if support >= 3 and sign >= 3 and rra < 0.10 and cv < 1.25:
            without_cct += 1
    except (ValueError, TypeError, KeyError):
        skipped_cct += 1

# Scenario: without logFC_CV
without_cv = 0
skipped_cv = 0
for r in tier_b_genes:
    try:
        support = int(r["support_n"])
        sign = int(r["sign_consistency_n"])
        rra = float(r["best_rra_fdr"])
        cct = float(r["best_cct_fdr"])
        if support >= 3 and sign >= 3 and rra < 0.10 and cct < 0.10:
            without_cv += 1
    except (ValueError, TypeError, KeyError):
        skipped_cv += 1

# Scenario: support ≥ 2 (relaxed)
support2 = 0
skipped_s2 = 0
for r in tier_b_genes:
    try:
        support = int(r["support_n"])
        sign = int(r["sign_consistency_n"])
        rra = float(r["best_rra_fdr"])
        cct = float(r["best_cct_fdr"])
        cv = float(r["logFC_CV"])
        if support >= 2 and sign >= 2 and rra < 0.10 and cct < 0.10 and cv < 1.25:
            support2 += 1
    except (ValueError, TypeError, KeyError):
        skipped_s2 += 1

current = len(tier_b_genes)
print(f"\n{'场景':<40} {'Tier B 基因数':>12} {'变化':>10}")
print(f"{'─'*65}")
print(f"{'当前 (所有 gates)':<40} {current:>12} {'—':>10}")
print(f"{'去掉 RRA gate (仅 CCT)':<40} {without_rra:>12} {without_rra-current:>+10}")
print(f"{'去掉 CCT gate (仅 RRA)':<40} {without_cct:>12} {without_cct-current:>+10}")
print(f"{'去掉 logFC_CV gate':<40} {without_cv:>12} {without_cv-current:>+10}")
print(f"{'放宽 support → 2':<40} {support2:>12} {support2-current:>+10}")
total_skipped = skipped_rra + skipped_cct + skipped_cv + skipped_s2
if total_skipped > 0:
    print(f"\n  注意: {total_skipped} 行因字段解析失败被跳过 (RRA={skipped_rra}, CCT={skipped_cct}, CV={skipped_cv}, S2={skipped_s2})")

# Also: how many current tier_B genes fail EACH specific gate?
fail_rra = sum(1 for r in tier_b_genes if float(r["best_rra_fdr"]) >= 0.10)
fail_cct = sum(1 for r in tier_b_genes if float(r["best_cct_fdr"]) >= 0.10)
fail_cv = sum(1 for r in tier_b_genes if float(r["logFC_CV"]) >= 1.25)
fail_support3 = sum(1 for r in tier_b_genes if int(r["support_n"]) < 3)
fail_sign3 = sum(1 for r in tier_b_genes if int(r["sign_consistency_n"]) < 3)

print(f"\n当前 Tier B 基因中，各 gate 单独失败数:")
print(f"  RRA FDR ≥ 0.10:          {fail_rra}")
print(f"  CCT FDR ≥ 0.10:          {fail_cct}")
print(f"  logFC_CV ≥ 1.25:         {fail_cv}")
print(f"  support_n < 3:           {fail_support3}")
print(f"  sign_consistency_n < 3:  {fail_sign3}")

# Which blocker b categories affect tier B genes?
print(f"\n当前 Tier B 基因的 tier_blocker_b 分布:")
blocker_b = Counter(r["tier_blocker_b"] for r in tier_b_genes)
for k, v in blocker_b.most_common(5):
    print(f"  {k}: {v}")

