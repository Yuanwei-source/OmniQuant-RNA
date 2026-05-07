#!/usr/bin/env Rscript
# workflow/scripts/benchmark_statistics.R
# Statistical tests + confidence intervals for all benchmark outputs.
# Read-only: loads existing files, computes stats, writes output.
# No downloads, no pipeline runs.
#
# Output: results/benchmark/statistical_tests.tsv

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

# ── Constants ──────────────────────────────────────────────────────────────────

DEG_PATH    <- "results/benchmark/annotation_degradation.tsv"
STAB_PATH   <- "results/benchmark/subsampling_stability.tsv"
ABL_PATH    <- "results/benchmark/ablation_summary.tsv"
OUT_PATH    <- "results/benchmark/statistical_tests.tsv"

# ── Helpers ────────────────────────────────────────────────────────────────────

cohens_d_paired <- function(x, y) {
  # Cohen's d for paired samples: mean(diff) / sd(diff)
  diff <- x - y
  mean_diff <- mean(diff, na.rm = TRUE)
  sd_diff <- sd(diff, na.rm = TRUE)
  if (is.na(sd_diff) || sd_diff == 0) return(NA_real_)
  mean_diff / sd_diff
}

std_error <- function(x) {
  sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
}

build_result_row <- function(comparison, metric, value_a, value_b,
                             effect_size, ci_lower, ci_upper, p_value) {
  tibble(
    comparison  = comparison,
    metric      = metric,
    value_a     = as.character(value_a),
    value_b     = as.character(value_b),
    effect_size = as.character(effect_size),
    ci_lower    = as.character(ci_lower),
    ci_upper    = as.character(ci_upper),
    p_value     = as.character(p_value)
  )
}

# ═══════════════════════════════════════════════════════════════════════════════
# A. DEGRADATION BENCHMARK — paired comparisons at 50% degradation level
# ═══════════════════════════════════════════════════════════════════════════════

cat("── A. Degradation benchmark ──\n")

deg <- read_tsv(DEG_PATH, show_col_types = FALSE)

# Focus on 50% degradation level
deg50 <- deg %>% filter(level == 0.50)

results_list <- list()
modes <- c("random_drop", "length_biased", "expression_biased", "transcript_level")
quantifiers <- c("featurecounts", "stringtie", "salmon", "kallisto")

for (m in modes) {
  dm <- deg50 %>% filter(mode == !!m)

  # Consensus data
  consensus_vals <- dm %>% filter(method == "consensus") %>%
    arrange(seed) %>% pull(global_recall)

  for (q in quantifiers) {
    q_vals <- dm %>% filter(method == !!q) %>%
      arrange(seed) %>% pull(global_recall)

    if (length(consensus_vals) < 2 || length(q_vals) < 2) {
      cat(sprintf("  WARNING: %s at %s level=0.5: insufficient data\n", m, q))
      next
    }

    # Align by seed (already sorted)
    paired <- complete.cases(consensus_vals, q_vals)
    if (sum(paired) < 2) {
      cat(sprintf("  WARNING: %s at %s level=0.5: no complete pairs\n", m, q))
      next
    }

    c_vals <- consensus_vals[paired]
    q_vals <- q_vals[paired]

    # Skip if data are essentially constant (sd of differences near zero)
    diff_sd <- sd(c_vals - q_vals, na.rm = TRUE)
    if (is.na(diff_sd) || diff_sd < 1e-10) {
      cat(sprintf("  %s vs %s at %s level=0.5: constant differences, skipping t-test\n",
                  "consensus", q, m))
      row <- build_result_row(
        comparison  = sprintf("degradation_%s_50pct", m),
        metric      = "global_recall",
        value_a     = sprintf("consensus=%.4f", mean(c_vals, na.rm = TRUE)),
        value_b     = sprintf("%s=%.4f", q, mean(q_vals, na.rm = TRUE)),
        effect_size = "— (constant diff)",
        ci_lower    = "—",
        ci_upper    = "—",
        p_value     = "—"
      )
      results_list[[length(results_list) + 1]] <- row
      next
    }

    # Paired t-test
    tt <- tryCatch(
      t.test(c_vals, q_vals, paired = TRUE),
      error = function(e) NULL
    )
    if (is.null(tt)) {
      cat(sprintf("  %s vs %s at %s: t.test failed, skipping\n", "consensus", q, m))
      next
    }

    # Cohen's d
    d <- cohens_d_paired(c_vals, q_vals)

    comp_name <- sprintf("degradation_%s_50pct", m)
    row <- build_result_row(
      comparison  = comp_name,
      metric      = "global_recall",
      value_a     = sprintf("consensus=%.4f", mean(c_vals, na.rm = TRUE)),
      value_b     = sprintf("%s=%.4f", q, mean(q_vals, na.rm = TRUE)),
      effect_size = sprintf("%.3f", d),
      ci_lower    = sprintf("%.4f", tt$conf.int[1]),
      ci_upper    = sprintf("%.4f", tt$conf.int[2]),
      p_value     = sprintf("%.4f", tt$p.value)
    )
    results_list[[length(results_list) + 1]] <- row
    cat(sprintf("  %s vs %s: d=%.3f CI=[%.4f,%.4f] p=%.4f\n",
                "consensus", q, d,
                tt$conf.int[1], tt$conf.int[2], tt$p.value))
  }
}

# Also compare featureCounts vs each other quantifier at 50%
for (m in modes) {
  dm <- deg50 %>% filter(mode == !!m)
  fc_vals <- dm %>% filter(method == "featurecounts") %>%
    arrange(seed) %>% pull(global_recall)

  for (q in setdiff(quantifiers, "featurecounts")) {
    q_vals <- dm %>% filter(method == !!q) %>%
      arrange(seed) %>% pull(global_recall)

    if (length(fc_vals) < 2 || length(q_vals) < 2) next
    paired <- complete.cases(fc_vals, q_vals)
    if (sum(paired) < 2) next

    f_vals <- fc_vals[paired]
    qq_vals <- q_vals[paired]

    # Skip if data are essentially constant (sd of differences near zero)
    diff_sd <- sd(f_vals - qq_vals, na.rm = TRUE)
    if (is.na(diff_sd) || diff_sd < 1e-10) {
      cat(sprintf("  %s vs %s at %s: constant data, skipping paired t-test\n",
                  "featurecounts", q, m))
      row <- build_result_row(
        comparison  = sprintf("degradation_%s_50pct", m),
        metric      = "global_recall",
        value_a     = sprintf("featurecounts=%.4f", mean(f_vals, na.rm = TRUE)),
        value_b     = sprintf("%s=%.4f", q, mean(qq_vals, na.rm = TRUE)),
        effect_size = "— (constant diff)",
        ci_lower    = "—",
        ci_upper    = "—",
        p_value     = "—"
      )
      results_list[[length(results_list) + 1]] <- row
      next
    }

    tt <- tryCatch(
      t.test(f_vals, qq_vals, paired = TRUE),
      error = function(e) NULL
    )
    if (is.null(tt)) {
      cat(sprintf("  %s vs %s at %s: t.test failed, skipping\n",
                  "featurecounts", q, m))
      next
    }

    d <- cohens_d_paired(f_vals, qq_vals)

    comp_name <- sprintf("degradation_%s_50pct", m)
    row <- build_result_row(
      comparison  = comp_name,
      metric      = "global_recall",
      value_a     = sprintf("featurecounts=%.4f", mean(f_vals, na.rm = TRUE)),
      value_b     = sprintf("%s=%.4f", q, mean(qq_vals, na.rm = TRUE)),
      effect_size = sprintf("%.3f", d),
      ci_lower    = sprintf("%.4f", tt$conf.int[1]),
      ci_upper    = sprintf("%.4f", tt$conf.int[2]),
      p_value     = sprintf("%.4f", tt$p.value)
    )
    results_list[[length(results_list) + 1]] <- row
    cat(sprintf("  %s vs %s: d=%.3f CI=[%.4f,%.4f] p=%.4f\n",
                "featurecounts", q, d,
                tt$conf.int[1], tt$conf.int[2], tt$p.value))
  }
}

# ═══════════════════════════════════════════════════════════════════════════════
# B. SUBSAMPLING — stability comparisons via Wilcoxon signed-rank
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n── B. Subsampling stability ──\n")

stab <- read_tsv(STAB_PATH, show_col_types = FALSE)

# Check that expected columns exist
expected_stab_cols <- c("stability_consensus", "stability_voting",
                        "stability_featurecounts", "stability_salmon",
                        "stability_stringtie", "stability_kallisto")
missing_cols <- setdiff(expected_stab_cols, colnames(stab))
if (length(missing_cols) > 0) {
  cat(sprintf("  WARNING: missing columns in subsampling_stability.tsv: %s\n",
              paste(missing_cols, collapse = ", ")))
}

consensus_stab <- stab$stability_consensus
voting_stab    <- stab$stability_voting

# Consensus vs each quantifier
stab_quantifiers <- c("featurecounts", "salmon", "stringtie", "kallisto")
stab_methods <- list(
  voting          = "stability_voting",
  featurecounts   = "stability_featurecounts",
  salmon          = "stability_salmon",
  stringtie       = "stability_stringtie",
  kallisto        = "stability_kallisto"
)

for (nm in names(stab_methods)) {
  col <- stab_methods[[nm]]
  if (!col %in% colnames(stab)) {
    cat(sprintf("  WARNING: column %s not found, skipping\n", col))
    next
  }

  vals <- stab[[col]]
  paired <- complete.cases(consensus_stab, vals)
  if (sum(paired) < 3) {
    cat(sprintf("  WARNING: insufficient paired data for %s\n", nm))
    next
  }

  c_val <- consensus_stab[paired]
  o_val <- vals[paired]

  # Wilcoxon signed-rank test (paired)
  wt <- wilcox.test(c_val, o_val, paired = TRUE, exact = FALSE,
                    conf.int = TRUE)

  # Effect size: median difference
  med_diff <- median(c_val - o_val, na.rm = TRUE)

  # For effect size, use r = Z / sqrt(N)
  # A non-parametric approximation of Cohen's d: use rank-biserial correlation
  n_pairs <- sum(paired)
  if (wt$statistic <= n_pairs * (n_pairs + 1) / 4) {
    r_effect <- 1 - (2 * wt$statistic) / (n_pairs * (n_pairs + 1))
  } else {
    # Use the complement
    v_star <- n_pairs * (n_pairs + 1) / 2 - wt$statistic
    r_effect <- (2 * v_star) / (n_pairs * (n_pairs + 1)) - 1
  }

  comp_name <- sprintf("stability_consensus_vs_%s", nm)
  row <- build_result_row(
    comparison  = comp_name,
    metric      = "gene_stability",
    value_a     = sprintf("consensus=%.4f", mean(c_val, na.rm = TRUE)),
    value_b     = sprintf("%s=%.4f", nm, mean(o_val, na.rm = TRUE)),
    effect_size = sprintf("%.3f (rank-biserial r)", r_effect),
    ci_lower    = sprintf("%.4f", wt$conf.int[1]),
    ci_upper    = sprintf("%.4f", wt$conf.int[2]),
    p_value     = sprintf("%.4f", wt$p.value)
  )
  results_list[[length(results_list) + 1]] <- row
  cat(sprintf("  consensus vs %s: median_diff=%.4f r=%.3f CI=[%.4f,%.4f] p=%.4f\n",
              nm, med_diff, r_effect,
              wt$conf.int[1], wt$conf.int[2], wt$p.value))
}

# Voting vs featureCounts (direct comparison)
if ("stability_voting" %in% colnames(stab) &&
    "stability_featurecounts" %in% colnames(stab)) {
  voting_vals <- stab$stability_voting
  fc_stab_vals <- stab$stability_featurecounts
  paired_vf <- complete.cases(voting_vals, fc_stab_vals)
  if (sum(paired_vf) >= 3) {
    v_val <- voting_vals[paired_vf]
    f_val <- fc_stab_vals[paired_vf]
    wt_vf <- wilcox.test(v_val, f_val, paired = TRUE, exact = FALSE,
                         conf.int = TRUE)
    comp_name <- "stability_voting_vs_featurecounts"
    row <- build_result_row(
      comparison  = comp_name,
      metric      = "gene_stability",
      value_a     = sprintf("voting=%.4f", mean(v_val, na.rm = TRUE)),
      value_b     = sprintf("featurecounts=%.4f", mean(f_val, na.rm = TRUE)),
      effect_size = "—",
      ci_lower    = sprintf("%.4f", wt_vf$conf.int[1]),
      ci_upper    = sprintf("%.4f", wt_vf$conf.int[2]),
      p_value     = sprintf("%.4f", wt_vf$p.value)
    )
    results_list[[length(results_list) + 1]] <- row
    cat(sprintf("  voting vs featurecounts: median_diff=%.4f CI=[%.4f,%.4f] p=%.4f\n",
                median(v_val - f_val, na.rm = TRUE),
                wt_vf$conf.int[1], wt_vf$conf.int[2], wt_vf$p.value))
  }
}

# ═══════════════════════════════════════════════════════════════════════════════
# C. ABLATION — percent changes from baseline
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n── C. Ablation study ──\n")

ablation <- read_tsv(ABL_PATH, show_col_types = FALSE)

# Extract values by config
baseline_row  <- ablation %>% filter(config == "baseline_single")
namespace_row <- ablation %>% filter(config == "+namespace")
consensus_row <- ablation %>% filter(config == "+consensus")
decontam_row  <- ablation %>% filter(config == "+decontam_full")

if (nrow(baseline_row) == 1 && nrow(consensus_row) == 1) {
  baseline_deg  <- baseline_row$n_deg
  consensus_deg <- consensus_row$n_deg
  pct_reduction <- (consensus_deg - baseline_deg) / baseline_deg * 100

  baseline_stab <- baseline_row$stability_mean
  consensus_stab <- consensus_row$stability_mean
  stab_change <- (consensus_stab - baseline_stab) / baseline_stab * 100

  row <- build_result_row(
    comparison  = "ablation_deg_reduction",
    metric      = "percent_change",
    value_a     = sprintf("baseline=%d", as.integer(baseline_deg)),
    value_b     = sprintf("consensus=%d", as.integer(consensus_deg)),
    effect_size = sprintf("%.1f%%", pct_reduction),
    ci_lower    = "—",
    ci_upper    = "—",
    p_value     = "—"
  )
  results_list[[length(results_list) + 1]] <- row
  cat(sprintf("  DEG reduction: baseline=%d → consensus=%d (%.1f%%)\n",
              as.integer(baseline_deg), as.integer(consensus_deg), pct_reduction))

  row2 <- build_result_row(
    comparison  = "ablation_stability_change",
    metric      = "percent_change",
    value_a     = sprintf("baseline=%.4f", baseline_stab),
    value_b     = sprintf("consensus=%.4f", consensus_stab),
    effect_size = sprintf("%.2f%%", stab_change),
    ci_lower    = "—",
    ci_upper    = "—",
    p_value     = "—"
  )
  results_list[[length(results_list) + 1]] <- row2
  cat(sprintf("  Stability change: baseline=%.4f → consensus=%.4f (%.2f%%)\n",
              baseline_stab, consensus_stab, stab_change))

  # Decontam host read percentage
  if (nrow(decontam_row) == 1) {
    host_pct <- decontam_row$host_reads_pct
    row3 <- build_result_row(
      comparison  = "ablation_host_reads",
      metric      = "percent",
      value_a     = "no_decontam=100%",
      value_b     = sprintf("decontam=%.2f%%", host_pct),
      effect_size = sprintf("%.1f%% host retained", host_pct),
      ci_lower    = "—",
      ci_upper    = "—",
      p_value     = "—"
    )
    results_list[[length(results_list) + 1]] <- row3
    cat(sprintf("  Host reads after decontam: %.2f%%\n", host_pct))
  }
} else {
  cat("  WARNING: could not find required rows in ablation_summary.tsv\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# D. WRITE OUTPUT
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n── Writing output ──\n")

dir.create(dirname(OUT_PATH), showWarnings = FALSE, recursive = TRUE)

if (length(results_list) > 0) {
  output <- bind_rows(results_list)
  write_tsv(output, OUT_PATH)
  cat(sprintf("  Wrote %d rows to %s\n", nrow(output), OUT_PATH))
} else {
  cat("  No results generated, writing empty output\n")
  output <- tibble(
    comparison  = character(),
    metric      = character(),
    value_a     = character(),
    value_b     = character(),
    effect_size = character(),
    ci_lower    = character(),
    ci_upper    = character(),
    p_value     = character()
  )
  write_tsv(output, OUT_PATH)
}

cat("\nStatistical tests complete.\n")
