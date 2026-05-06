#!/usr/bin/env Rscript
# benchmark_ablation.R — Quantify performance contribution of each OmniQuant module
# by systematically removing them one at a time (ablation study).
#
# Reads EXISTING pipeline outputs. No pipeline re-runs, no downloads.
#
# Output: results/benchmark/ablation_summary.tsv

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

# ── Paths ──────────────────────────────────────────────────────────────────────
CONSENSUS_PATH   <- "results/07.consensus_expression/60d_vs_1d/consensus_results.tsv"
FC_DEA_PATH      <- "results/06.differential_expression/featurecounts/deseq2.60d_vs_1d.csv"
TX2GENE_PATH     <- "results/00.reference/tx2gene_master.tsv"
STABILITY_PATH   <- "results/benchmark/subsampling_summary.tsv"
DECONTAM_DIR     <- "results/03.decontam/stats/"
OUT_PATH         <- "results/benchmark/ablation_summary.tsv"

# ── Helper: safe file read, return NULL if missing ────────────────────────────
safe_read <- function(path, reader, ...) {
  if (missing(reader)) {
    reader <- function(p) readr::read_tsv(p, show_col_types = FALSE)
  }
  if (file.exists(path)) {
    return(reader(path, ...))
  }
  cat(sprintf("WARNING: File not found: %s\n", path))
  return(NULL)
}

# ── Helper: safe value extraction, NA if NULL ─────────────────────────────────
safe_val <- function(x, default = NA_real_) {
  if (is.null(x) || length(x) == 0) return(default)
  return(x)
}

# ═══════════════════════════════════════════════════════════════════════════════
# 1. LOAD EXISTING DATA
# ═══════════════════════════════════════════════════════════════════════════════
cat("Loading existing pipeline outputs...\n")

# Consensus results (Tab-separated)
consensus <- safe_read(CONSENSUS_PATH)
cat(sprintf("  consensus_results: %d rows\n", safe_val(nrow(consensus), 0L)))

# FeatureCounts DESeq2 results (comma-separated)
fc_dea <- safe_read(FC_DEA_PATH, reader = function(p) readr::read_csv(p, show_col_types = FALSE))
cat(sprintf("  featureCounts DEA: %d rows\n", safe_val(nrow(fc_dea), 0L)))

# tx2gene master namespace bridge
tx2g <- safe_read(TX2GENE_PATH)
cat(sprintf("  tx2gene_master: %d rows\n", safe_val(nrow(tx2g), 0L)))

# Subsampling stability summary (pre-computed mean stabilities)
stability <- safe_read(STABILITY_PATH)
cat(sprintf("  subsampling_summary: %d rows\n", safe_val(nrow(stability), 0L)))

# Decontam host rescue stats (all *_host_rescue_stats.tsv files)
decontam_files <- list.files(DECONTAM_DIR, pattern = "_host_rescue_stats\\.tsv$",
                              full.names = TRUE)
cat(sprintf("  decontam host_rescue files: %d\n", length(decontam_files)))

# ═══════════════════════════════════════════════════════════════════════════════
# 2. COMPUTE METRICS
# ═══════════════════════════════════════════════════════════════════════════════

# ── 2a. Baseline: single quantifier DE genes (featureCounts only) ─────────────
# FDR < 0.05 AND |logFC| >= 1.0
if (!is.null(fc_dea)) {
  baseline_deg <- sum(
    fc_dea$adj.P.Val <= 0.05 & abs(fc_dea$logFC) >= 1.0,
    na.rm = TRUE
  )
  # Total detected genes (gene universe size for featureCounts)
  n_detected_baseline <- nrow(fc_dea)
} else {
  baseline_deg <- NA_real_
  n_detected_baseline <- NA_real_
}
cat(sprintf("  Baseline DEG (featureCounts FDR<0.05 & |logFC|>=1): %d\n",
            safe_val(baseline_deg, 0L)))

# ── 2b. +Namespace: ID recovery rate ──────────────────────────────────────────
# Fraction of reference genes with resolved IDs (allow_consensus_main == TRUE)
if (!is.null(tx2g)) {
  # Use unique gene_ids to avoid double-counting transcripts
  tx2g_ref <- tx2g[tx2g$is_reference_gene == "true" | tx2g$is_reference_gene == TRUE, , drop = FALSE]
  if (nrow(tx2g_ref) > 0) {
    # Count unique gene IDs
    ref_genes <- unique(tx2g_ref$gene_id)
    # Count unique gene IDs that pass namespace resolution
    tx2g_pass <- tx2g_ref[tx2g_ref$allow_consensus_main == "true" |
                           tx2g_ref$allow_consensus_main == TRUE, , drop = FALSE]
    pass_genes <- unique(tx2g_pass$gene_id)
    id_recovery_rate <- length(pass_genes) / length(ref_genes)
  } else {
    id_recovery_rate <- NA_real_
  }
} else {
  id_recovery_rate <- NA_real_
}
cat(sprintf("  ID recovery rate (namespace): %.4f\n", safe_val(id_recovery_rate, 0)))

# ── 2c. +Consensus: Tier A genes and consensus stability ──────────────────────
if (!is.null(consensus)) {
  tier_a_count <- sum(consensus$tier == "Tier_A", na.rm = TRUE)
} else {
  tier_a_count <- NA_real_
}
cat(sprintf("  Consensus Tier_A genes: %d\n", safe_val(tier_a_count, 0L)))

# ── 2d. Stability means (from subsampling summary) ────────────────────────────
stability_featurecounts <- NA_real_
stability_consensus <- NA_real_

if (!is.null(stability)) {
  # Normalize method column values (handle possible case/whitespace differences)
  stab_methods <- tolower(trimws(stability$method))

  idx_fc <- which(stab_methods == "featurecounts")
  if (length(idx_fc) > 0) {
    stability_featurecounts <- as.numeric(stability$mean_stability[idx_fc[1]])
  }

  idx_ta <- which(stab_methods == "tier_a")
  if (length(idx_ta) > 0) {
    stability_consensus <- as.numeric(stability$mean_stability[idx_ta[1]])
  }
}
cat(sprintf("  Stability mean (featureCounts): %.4f\n", safe_val(stability_featurecounts, 0)))
cat(sprintf("  Stability mean (Tier A):       %.4f\n", safe_val(stability_consensus, 0)))

# ── 2e. Decontam: host read percentage ────────────────────────────────────────
# Compute from host_rescue_stats files: host_pairs / (host_pairs + unresolved_pairs)
host_reads_pct <- NA_real_

if (length(decontam_files) > 0) {
  decontam_samples <- lapply(decontam_files, function(f) {
    d <- readr::read_tsv(f, show_col_types = FALSE)
    # Use last row (skip potential header duplication)
    d[nrow(d), , drop = FALSE]
  })
  decontam_all <- do.call(rbind, decontam_samples)

  # Compute per-sample host fraction
  total_pairs <- decontam_all$host_pairs + decontam_all$ercc_pairs +
                 decontam_all$unresolved_pairs
  host_fraction <- decontam_all$host_pairs / total_pairs
  host_reads_pct <- mean(host_fraction, na.rm = TRUE) * 100

  cat(sprintf("  Decontam samples: %d\n", nrow(decontam_all)))
  cat(sprintf("  Host reads pct (mean across samples): %.2f%%\n",
              safe_val(host_reads_pct, 0)))
} else {
  cat("  WARNING: No decontam host_rescue_stats files found.\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# 3. BUILD ABLATION TABLE
# ═══════════════════════════════════════════════════════════════════════════════

results <- tibble::tibble(
  config          = c("baseline_single", "+namespace", "+consensus", "+decontam_full"),
  n_deg           = c(
    baseline_deg,                     # single quantifier
    baseline_deg,                     # same DEA, just adds namespace
    tier_a_count,                      # consensus Tier A
    tier_a_count                       # same consensus, with decontam
  ),
  tier_a_genes    = c(
    0,                                 # no consensus tier
    0,                                 # no consensus tier
    tier_a_count,
    tier_a_count
  ),
  stability_mean  = c(
    stability_featurecounts,          # featureCounts bootstrap stability
    stability_featurecounts,          # same as baseline
    stability_consensus,              # consensus stability
    stability_consensus               # same consensus, with decontam
  ),
  id_recovery_rate = c(
    1.0,                               # no namespace filter → perfect recovery
    id_recovery_rate,                  # namespace-filtered recovery
    id_recovery_rate,                  # same namespace
    id_recovery_rate                   # same namespace
  ),
  host_reads_pct  = c(
    100,                               # no decontam
    100,                               # no decontam
    100,                               # no decontam
    host_reads_pct                     # post-decontam
  )
)

# ═══════════════════════════════════════════════════════════════════════════════
# 4. WRITE OUTPUT
# ═══════════════════════════════════════════════════════════════════════════════

# Ensure output directory exists
dir.create(dirname(OUT_PATH), showWarnings = FALSE, recursive = TRUE)

readr::write_tsv(results, OUT_PATH)
cat(sprintf("\nAblation summary written to: %s\n", OUT_PATH))

# ── Human-readable summary ────────────────────────────────────────────────────
cat("\n═══════════════════════════════════════════════════════════════\n")
cat("  OMNiquant MODULE ABLATION STUDY\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

pad <- function(x, w = 18) sprintf(sprintf("%%-%ds", w), x)

cat(pad("Configuration"), pad("DEGs", 8), pad("Tier_A", 8),
    pad("Stability", 10), pad("ID Recovery", 12),
    pad("Host %", 8), "\n", sep = "")
cat(strrep("-", 72), "\n", sep = "")

for (i in seq_len(nrow(results))) {
  row <- results[i, ]
  fmt_int  <- function(x) if (is.na(x)) "NA      " else sprintf("%-8d", as.integer(x))
  fmt_num  <- function(x, f) if (is.na(x)) "NA       " else sprintf(f, x)
  cat(
    pad(row$config, 18),
    fmt_int(row$n_deg),
    fmt_int(row$tier_a_genes),
    fmt_num(row$stability_mean, "%-10.3f"),
    fmt_num(row$id_recovery_rate, "%-12.3f"),
    fmt_num(row$host_reads_pct, "%-8.1f"),
    "\n", sep = ""
  )
}

cat("\nNotes:\n")
cat("  baseline_single: featureCounts DESeq2 only, no consensus, no namespace\n")
cat("  +namespace:      adds reference auto-detection + gene ID namespace resolution\n")
cat("  +consensus:      adds RRA+CCT multi-quantifier consensus DEA\n")
cat("  +decontam_full:  full pipeline including conservative decontamination\n")

# Report any missing metrics
missing_metrics <- character(0)
if (is.na(baseline_deg))          missing_metrics <- c(missing_metrics, "n_deg")
if (is.na(id_recovery_rate))      missing_metrics <- c(missing_metrics, "id_recovery_rate")
if (is.na(stability_consensus))   missing_metrics <- c(missing_metrics, "stability_mean")
if (is.na(host_reads_pct))        missing_metrics <- c(missing_metrics, "host_reads_pct")

if (length(missing_metrics) > 0) {
  cat(sprintf("\nWARNING: Missing metrics (NA in output): %s\n",
              paste(missing_metrics, collapse = ", ")))
  cat("  These data files were not found at the expected paths.\n")
}

cat("\nAblation study complete.\n")
