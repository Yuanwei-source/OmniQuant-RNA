# Learnings: Drosophila Subsampling Benchmark

## Count matrix formats
- featureCounts: column `Geneid` with `gene:FBgn...` IDs, integer counts
- Salmon: column `gene_id` with `gene:FBgn...` IDs, float counts (tximport)
- StringTie: column `Gene_ID` with `gene:FBgn...` IDs, integer counts
- Kallisto: column `target_id` with `transcript:FBtr...` IDs, plus `length`/`eff_length` metadata columns, float counts

## tx2gene_master
- Column `allow_consensus_main` has lowercase `true`/`false` — R's `as.logical()` returns NA
- Must use `as_flag()` helper (tolower + %in% c("true","1","yes"))

## Consensus functions replicated
- `choose_direction()`: 2-arg version (up_support, down_support) — no p-value tiebreaker
- `assign_tier()`: support_n < 2 → "unclassified"
- `compute_rra_scores()`: exact=FALSE (Bonferroni-corrected)
- CCT thresholds set to 1.0 per plan → RRA-driven consensus

## Samples
- 3 reps per group (60d and 1d), 3 choose 2 = 3 combos per group, 3×3 = 9 subsets
- DESeq2: contrast=c("group", "60d", "1d") — 60d is treatment (numerator)

## Edge cases handled
- Genes with all-zero counts in subset → removed before DESeq2
- DESeq2 failures → tryCatch returns empty tibble
- Salmon float counts → rounded to integer
- replace_na not available in dplyr alone → used coalesce()

## Annotation degradation benchmark (2026-05-06)

### GFF namespace mismatch
- GFF3: `gene_id=FBgnXXXX` (no prefix), `ID=gene:FBgnXXXX` (with prefix)
- Consensus/featureCounts: use `gene:FBgnXXXX` format
- Must extract `ID` attribute for gene entries to match consensus namespace
- Standardize: prepend "gene:" to gene_id if not already prefixed

### Package loading
- `box` package not always available → added fallback to `library()` calls
- dplyr `group_by()`/`summarise()`/`n()` — use `dplyr::` prefix for portability
- `dp$group_by()` inside `%>%` works with box::use but not with library()

### dplyr column manipulation pitfalls
- `dea_df[mask, id_col]` returns data.frame (list) when used with dplyr tibbles
- Must use `dea_df[[id_col]][mask]` for vector extraction
- Column assignment with string names: `dea_df[mask, "logFC"] <- NA_real_` works fine

### Degradation simulation insights
- Transcript-level at 75%: 0 genes lost all isoforms (Drosophila has ~3 isoforms/gene avg)
- Expression-biased drop: highest global_recall (DE genes tend to be highly expressed)
- ID corruption: precision drops proportionally (corrupted IDs create false detections)
- Random drop vs length-biased: similar patterns (annotation-dependent quantifiers equally affected)

### Level normalization
- CLI accepts both fractions (0.75) and percentages (75)
- Auto-normalize: if any value > 1, divide by 100
- Levels displayed as `level * 100`% in output

### Output format verified
- Per-seed TSV: mode, level, method, detectable_recall, global_recall, precision, n_detected, n_reference, seed
- Summary TSV: mean ± SD across seeds
- 4 figures: global_recall, detectable_recall, precision, combined multi-panel

## Module ablation benchmark (2026-05-06)

### Data sources used
- `results/07.consensus_expression/60d_vs_1d/consensus_results.tsv` — 18,315 genes, 4,740 Tier_A
- `results/06.differential_expression/featurecounts/deseq2.60d_vs_1d.csv` — comma-separated (CSV), not TSV!
- `results/00.reference/tx2gene_master.tsv` — 71,320 rows (35,660 reference + 35,660 stringtie)
- `results/benchmark/subsampling_summary.tsv` — pre-computed mean_stability per method
- `results/03.decontam/stats/*_host_rescue_stats.tsv` — per-sample host/unresolved pairs

### Key data format differences
- Consensus file: tab-separated (TSV)
- FeatureCounts DEA file: comma-separated (CSV) — needs `readr::read_csv()` not `read_tsv()`
- Must handle reader type per file in safe_read wrapper

### Decontam stats interpretation
- `host_rescue_stats.tsv` has `host_pairs`, `ercc_pairs`, `unresolved_pairs`
- All samples: ercc_pairs=0, classified_pairs=0, nontarget=0, uncertain=0
- host_reads_pct = mean(host_pairs / (host_pairs + unresolved_pairs)) * 100 → ~65%
- `project_decontam_summary.tsv` has `Rescued_ERCC` matching `host_pairs` — likely misnamed column

### ID recovery rate
- All 18,315 reference genes have `allow_consensus_main = TRUE` → id_recovery_rate = 1.0
- Gene-level counting (unique gene_ids) avoids double-counting from multiple transcripts
- `is_reference_gene` values: lowercase "true" string, not boolean

### Metrics computed
- n_deg (baseline): featureCounts DESeq2 FDR<0.05 & |logFC|>=1 → 5,995 genes
- tier_a_genes: consensus Tier_A count → 4,740 genes
- stability_mean: from subsampling_summary.tsv mean_stability column
  - featurecounts: 0.7856
  - tier_a: 0.7824
- id_recovery_rate: 1.0 (all ref genes pass namespace)
- host_reads_pct: 65.23% (mean across 6 samples from host_rescue_stats)

### Output format
- `results/benchmark/ablation_summary.tsv` — 6 columns: config, n_deg, tier_a_genes, stability_mean, id_recovery_rate, host_reads_pct
- 4 configurations: baseline_single, +namespace, +consensus, +decontam_full
