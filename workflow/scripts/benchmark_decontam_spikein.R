#!/usr/bin/env Rscript
# workflow/scripts/benchmark_decontam_spikein.R
# Decontam module benchmark — condition-specific microbial spike-in evaluation
# Four modes:
#   assess  — evaluate current contamination levels (→ decontam_assessment.txt)
#   spike   — create condition-specific spiked FASTQ
#   compare — analyze decontam ON/OFF DEG + decontam efficiency results
#   all     — run assess, then conditionally spike or compare

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

# ── Parse CLI args ──────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
opts <- list(
  mode            = "assess",
  microbe_fastq   = NULL,
  spike_rate_a    = 0.15,
  spike_rate_b    = 0.01,
  stat_dir        = "results/03.decontam/stats",
  samples_tsv     = "data/fastq/samples.tsv",
  clean_dir       = "results/03.decontam/clean",
  output_dir      = "results/benchmark",
  on_dir          = "results.backup.decontam_on",
  off_dir         = "results.backup.decontam_off",
  on_stats_dir    = NULL,
  off_stats_dir   = NULL,
  contrast        = NULL
)
i <- 1
while (i <= length(args)) {
  if (args[i] == "--mode")                { opts$mode <- args[i+1]; i <- i+2 }
  else if (args[i] == "--microbe-fastq")  { opts$microbe_fastq <- args[i+1]; i <- i+2 }
  else if (args[i] == "--spike-rate-a")   { opts$spike_rate_a <- as.numeric(args[i+1]); i <- i+2 }
  else if (args[i] == "--spike-rate-b")   { opts$spike_rate_b <- as.numeric(args[i+1]); i <- i+2 }
  else if (args[i] == "--stat-dir")       { opts$stat_dir <- args[i+1]; i <- i+2 }
  else if (args[i] == "--samples-tsv")    { opts$samples_tsv <- args[i+1]; i <- i+2 }
  else if (args[i] == "--clean-dir")      { opts$clean_dir <- args[i+1]; i <- i+2 }
  else if (args[i] == "--output-dir")     { opts$output_dir <- args[i+1]; i <- i+2 }
  else if (args[i] == "--on-dir")         { opts$on_dir <- args[i+1]; i <- i+2 }
  else if (args[i] == "--off-dir")        { opts$off_dir <- args[i+1]; i <- i+2 }
  else if (args[i] == "--on-stats-dir")   { opts$on_stats_dir <- args[i+1]; i <- i+2 }
  else if (args[i] == "--off-stats-dir")  { opts$off_stats_dir <- args[i+1]; i <- i+2 }
  else if (args[i] == "--contrast")       { opts$contrast <- args[i+1]; i <- i+2 }
  else { i <- i+1 }
}

dir.create(opts$output_dir, recursive = TRUE, showWarnings = FALSE)

# ── Shared helpers ──────────────────────────────────────────────────────────────

safe_read_tsv <- function(path, ...) {
  if (file.exists(path)) {
    tryCatch(read_tsv(path, show_col_types = FALSE, comment = "#", ...),
             error = function(e) { cat(sprintf("WARNING: could not read %s: %s\n", path, e$message)); NULL })
  } else {
    cat(sprintf("WARNING: file not found: %s\n", path))
    NULL
  }
}

collect_decontam_stats <- function(stat_dir) {
  # Returns per-sample fate breakdown from decontam decision_summary files
  if (!dir.exists(stat_dir)) {
    cat(sprintf("WARNING: decontam stats directory not found: %s\n", stat_dir))
    return(NULL)
  }
  decision_files <- list.files(stat_dir, pattern = "_decision_summary.tsv", full.names = TRUE)
  if (length(decision_files) == 0) {
    cat(sprintf("WARNING: no decision_summary files in %s\n", stat_dir))
    return(NULL)
  }
  rows <- lapply(decision_files, function(f) {
    d <- safe_read_tsv(f)
    if (is.null(d)) return(NULL)
    sample_id <- sub("_decision_summary.tsv$", "", basename(f))
    d$sample <- sample_id
    total <- sum(d$Rescued_Host, d$Rescued_ERCC, d$Flagged_Uncertain,
                 d$Removed_Tech, d$Removed_NonTarget, na.rm = TRUE)
    d %>% mutate(total_pairs = total,
                 host_pct     = round(100 * Rescued_Host / total, 2),
                 uncertain_pct = round(100 * Flagged_Uncertain / total, 2),
                 removed_pct  = round(100 * (Removed_Tech + Removed_NonTarget) / total, 2))
  })
  bind_rows(rows)
}

collect_classification_stats <- function(stat_dir) {
  if (!dir.exists(stat_dir)) return(NULL)
  class_files <- list.files(stat_dir, pattern = "_classification_stats.tsv", full.names = TRUE)
  if (length(class_files) == 0) return(NULL)
  rows <- lapply(class_files, function(f) {
    d <- safe_read_tsv(f)
    if (is.null(d)) return(NULL)
    d$sample <- sub("_classification_stats.tsv$", "", basename(f))
    d
  })
  bind_rows(rows)
}

# ══════════════════════════════════════════════════════════════════════════════════
# MODE 1: ASSESS — evaluate current contamination levels
# ══════════════════════════════════════════════════════════════════════════════════
if (opts$mode == "assess") {
  cat("\n=== Decontam Spike-in Benchmark: ASSESS MODE ===\n\n")

  stat_dir <- opts$stat_dir
  if (!dir.exists(stat_dir)) {
    cat("ERROR: Decontam stats not found at", stat_dir, "\n")
    cat("Run the full pipeline with decontam.enabled=true first.\n")
    quit(status = 1)
  }

  classification_files <- list.files(stat_dir, pattern = "_classification_stats.tsv", full.names = TRUE)
  rescue_files         <- list.files(stat_dir, pattern = "_host_rescue_stats.tsv", full.names = TRUE)
  decision_files       <- list.files(stat_dir, pattern = "_decision_summary.tsv", full.names = TRUE)

  total_host    <- 0
  total_nonhost <- 0
  total_classified   <- 0
  total_unclassified <- 0
  total_uncertain    <- 0
  per_sample <- list()

  for (f in classification_files) {
    d <- tryCatch(read_tsv(f, show_col_types = FALSE, comment = "#"), error = function(e) NULL)
    if (is.null(d)) next
    needed <- c("unclassified_pairs")
    if (!all(needed %in% names(d))) next
    s <- sub("_classification_stats.tsv$", "", basename(f))
    nonhost <- sum(d$unclassified_pairs, na.rm = TRUE)
    total_nonhost <- total_nonhost + nonhost
    per_sample[[s]] <- c(per_sample[[s]], nonhost_pairs = nonhost)
    if ("classified_pairs" %in% names(d))
      total_classified <- total_classified + sum(d$classified_pairs, na.rm = TRUE)
    if ("uncertain_pairs" %in% names(d))
      total_uncertain <- total_uncertain + sum(d$uncertain_pairs, na.rm = TRUE)
  }

  for (f in rescue_files) {
    d <- tryCatch(read_tsv(f, show_col_types = FALSE, comment = "#"), error = function(e) NULL)
    if (is.null(d)) next
    if (!"host_pairs" %in% names(d)) next
    s <- sub("_host_rescue_stats.tsv$", "", basename(f))
    host <- sum(d$host_pairs, na.rm = TRUE)
    total_host <- total_host + host
    per_sample[[s]] <- c(per_sample[[s]], host_pairs = host)
  }

  total_pairs  <- total_host + total_nonhost
  nonhost_pct  <- if (total_pairs > 0) round(100 * total_nonhost / total_pairs, 2) else 0
  host_pct     <- if (total_pairs > 0) round(100 * total_host / total_pairs, 2) else 0

  cat(sprintf("Total read pairs:    %d\n", total_pairs))
  cat(sprintf("Host read pairs:     %d (%.2f%%)\n", total_host, host_pct))
  cat(sprintf("Non-host pairs:      %d (%.2f%%)\n", total_nonhost, nonhost_pct))
  if (total_classified > 0)
    cat(sprintf("  Classified:        %d\n", total_classified))
  if (total_uncertain > 0)
    cat(sprintf("  Uncertain:         %d\n", total_uncertain))
  cat("\n")

  # Per-sample detail
  if (length(per_sample) > 0) {
    ps_df <- do.call(rbind, lapply(names(per_sample), function(s) {
      vals <- per_sample[[s]]
      data.frame(sample = s, host_pairs = vals["host_pairs"], nonhost_pairs = vals["nonhost_pairs"],
                 stringsAsFactors = FALSE)
    })) %>%
      mutate(total = host_pairs + nonhost_pairs,
             nonhost_pct = round(100 * nonhost_pairs / total, 2))
    cat("Per-sample breakdown:\n")
    print(ps_df, row.names = FALSE)
    write_tsv(ps_df, file.path(opts$output_dir, "decontam_per_sample_assessment.tsv"))
  }

  if (nonhost_pct > 5) {
    cat(sprintf("\nRECOMMENDATION: Non-host reads %.2f%% > 5%%. Use real data for decontam benchmark.\n", nonhost_pct))
    cat("  1. Run pipeline with decontam.enabled=true  → backup to results.backup.decontam_on\n")
    cat("  2. Run pipeline with decontam.enabled=false → backup to results.backup.decontam_off\n")
    cat("  3. Compare:\n")
    cat(sprintf("     Rscript workflow/scripts/benchmark_decontam_spikein.R --mode compare \\\n"))
    cat(sprintf("       --on-dir results.backup.decontam_on --off-dir results.backup.decontam_off \\\n"))
    cat(sprintf("       --contrast 60d_vs_1d --output-dir %s\n", opts$output_dir))
  } else {
    cat(sprintf("\nRECOMMENDATION: Non-host reads %.2f%% < 5%%. Spike-in needed.\n", nonhost_pct))
    cat("\nTo download Wolbachia reads for spike-in:\n")
    cat("  conda activate kingfisher\n")
    cat("  kingfisher get -r SRR1553797 -m ena-ftp -o results/benchmark/\n")
    cat("\nThen run spike-in:\n")
    cat(sprintf("  Rscript workflow/scripts/benchmark_decontam_spikein.R --mode spike \\\n"))
    cat("    --microbe-fastq results/benchmark/wolbachia.fastq.gz \\\n")
    cat("    --spike-rate-a 0.15 --spike-rate-b 0.01\n")
  }

  writeLines(c(
    sprintf("host_pairs=%d", total_host),
    sprintf("nonhost_pairs=%d", total_nonhost),
    sprintf("nonhost_pct=%.2f", nonhost_pct),
    sprintf("recommendation=%s", if (nonhost_pct > 5) "use_real_data" else "spike_in_needed")
  ), file.path(opts$output_dir, "decontam_assessment.txt"))

  cat(sprintf("\nAssessment saved to %s/decontam_assessment.txt\n", opts$output_dir))
}

# ══════════════════════════════════════════════════════════════════════════════════
# MODE 2: SPIKE — create condition-specific spiked FASTQ
# ══════════════════════════════════════════════════════════════════════════════════
if (opts$mode == "spike") {
  cat("\n=== Decontam Spike-in Benchmark: SPIKE MODE ===\n\n")

  if (is.null(opts$microbe_fastq) || !file.exists(opts$microbe_fastq)) {
    cat("ERROR: --microbe-fastq not specified or file not found\n")
    quit(status = 1)
  }

  if (!dir.exists(opts$clean_dir)) {
    cat(sprintf("ERROR: clean FASTQ directory not found: %s\n", opts$clean_dir))
    cat("Run the pipeline with decontam.enabled=true first to generate clean reads.\n")
    quit(status = 1)
  }

  samples <- safe_read_tsv(opts$samples_tsv)
  if (is.null(samples)) {
    cat(sprintf("ERROR: samples.tsv not found: %s\n", opts$samples_tsv))
    quit(status = 1)
  }
  if (!"group" %in% names(samples)) {
    colnames(samples)[1:4] <- c("sample", "fq1", "fq2", "group")
  }

  spiked_dir <- file.path(opts$output_dir, "spiked_fastq")
  dir.create(spiked_dir, recursive = TRUE, showWarnings = FALSE)

  cat(sprintf("Spike-in rates: Condition A=%.0f%%, Condition B=%.0f%%\n",
              opts$spike_rate_a * 100, opts$spike_rate_b * 100))
  cat(sprintf("Microbial source: %s\n", opts$microbe_fastq))
  cat("\n")

  # Check for seqtk
  seqtk_path <- Sys.which("seqtk")
  if (seqtk_path == "") {
    cat("ERROR: seqtk not found in PATH. Install with: conda install -c bioconda seqtk\n")
    quit(status = 1)
  }

  # Determine condition groups
  groups <- unique(samples$group)
  if (length(groups) < 2) {
    cat("ERROR: samples.tsv must have at least 2 groups for condition-specific spike-in\n")
    quit(status = 1)
  }
  group_a <- groups[1]
  group_b <- groups[2]
  cat(sprintf("Group assignments: '%s' → %.0f%% spike-in, '%s' → %.0f%% spike-in\n",
              group_a, opts$spike_rate_a * 100, group_b, opts$spike_rate_b * 100))

  # Count total microbial reads needed
  microbe_lines <- as.integer(system(sprintf("zcat -f %s | wc -l", shQuote(opts$microbe_fastq)),
                                    intern = TRUE))
  microbe_reads <- microbe_lines / 4
  cat(sprintf("Available microbial reads: %d\n", microbe_reads))

  new_samples <- list()
  genome_size <- 1e8  # rough estimate for seeding

  for (idx in seq_len(nrow(samples))) {
    s     <- samples$sample[idx]
    group <- samples$group[idx]
    rate  <- if (group == group_a) opts$spike_rate_a else opts$spike_rate_b

    r1_clean <- file.path(opts$clean_dir, paste0(s, "_R1_clean.fastq.gz"))
    r2_clean <- file.path(opts$clean_dir, paste0(s, "_R2_clean.fastq.gz"))

    if (!file.exists(r1_clean)) {
      cat(sprintf("WARNING: clean FASTQ not found for sample %s, skipping\n", s))
      next
    }

    host_lines <- as.integer(system(sprintf("zcat -f %s | wc -l", shQuote(r1_clean)), intern = TRUE))
    host_reads <- host_lines / 4
    spike_reads_needed <- round(host_reads * rate)

    if (spike_reads_needed > microbe_reads) {
      cat(sprintf("WARNING: sample %s needs %d spike reads but only %d available, capping\n",
                  s, spike_reads_needed, microbe_reads))
      spike_reads_needed <- microbe_reads
    }

    cat(sprintf("Sample %s (%s): host=%d reads, spike-in=%d reads (%.1f%%)\n",
                s, group, host_reads, spike_reads_needed, rate * 100))

    # Subsample microbial reads with seqtk
    seed <- as.integer(Sys.time()) %% .Machine$integer.max + idx
    microbe_r1_sub <- file.path(spiked_dir, paste0(s, "_microbe_R1.fastq.gz"))
    microbe_r2_sub <- file.path(spiked_dir, paste0(s, "_microbe_R2.fastq.gz"))

    # seqtk sample uses seed-based deterministic sampling; need paired-end mode
    cmd_sample <- sprintf(
      "seqtk sample -s%d %s %d 2>/dev/null | gzip -c > %s",
      seed, shQuote(opts$microbe_fastq), spike_reads_needed, shQuote(microbe_r1_sub)
    )
    ret <- system(cmd_sample)
    if (ret != 0) {
      cat(sprintf("ERROR: seqtk sample failed for sample %s\n", s))
      next
    }

    # If microbe FASTQ is paired, we need to handle pairs. For now, single-end spike-in:
    # duplicate as R2. If user provides paired microbe FASTQ, this should be adjusted.
    system(sprintf("cp %s %s", shQuote(microbe_r1_sub), shQuote(microbe_r2_sub)))

    # Concatenate host + spike
    r1_out <- file.path(spiked_dir, paste0(s, "_R1_spiked.fastq.gz"))
    r2_out <- file.path(spiked_dir, paste0(s, "_R2_spiked.fastq.gz"))

    system(sprintf("cat %s %s > %s", shQuote(r1_clean), shQuote(microbe_r1_sub), shQuote(r1_out)))
    system(sprintf("cat %s %s > %s", shQuote(r2_clean), shQuote(microbe_r2_sub), shQuote(r2_out)))

    new_samples[[length(new_samples) + 1]] <- data.frame(
      sample = s, fq1 = r1_out, fq2 = r2_out, group = group, stringsAsFactors = FALSE
    )
  }

  if (length(new_samples) == 0) {
    cat("ERROR: no samples were processed\n")
    quit(status = 1)
  }

  new_samples_df <- bind_rows(new_samples)
  new_samples_path <- file.path(opts$output_dir, "samples_spiked.tsv")
  write_tsv(new_samples_df, new_samples_path)

  cat(sprintf("\nCreated %d spiked FASTQ files in %s\n", nrow(new_samples_df) * 2, spiked_dir))
  cat(sprintf("Updated sample table: %s\n", new_samples_path))
  cat("\nNext steps:\n")
  cat("  1. cp config/config.yaml config/config_backup.yaml\n")
  cat("  2. Edit config.yaml: samples → new samples_spiked.tsv path\n")
  cat("  3. Run decontam OFF:  config.yaml → decontam.enabled=false → backup results\n")
  cat("  4. Run decontam ON:   config.yaml → decontam.enabled=true  → backup results\n")
  cat(sprintf("  5. Compare: Rscript workflow/scripts/benchmark_decontam_spikein.R --mode compare \\\n"))
  cat(sprintf("       --on-dir <on_results> --off-dir <off_results> --contrast <contrast> \\\n"))
  cat(sprintf("       --output-dir %s\n", opts$output_dir))
}

# ══════════════════════════════════════════════════════════════════════════════════
# MODE 3: COMPARE — analyze decontam ON/OFF consensus DEA + decontam efficiency
# ══════════════════════════════════════════════════════════════════════════════════
if (opts$mode == "compare") {
  cat("\n=== Decontam Spike-in Benchmark: COMPARE MODE ===\n\n")

  # Resolve contrast: if not specified, auto-detect from ON directory
  consensus_dir_on  <- file.path(opts$on_dir, "07.consensus_expression")
  consensus_dir_off <- file.path(opts$off_dir, "07.consensus_expression")

  if (!dir.exists(consensus_dir_on)) {
    cat(sprintf("ERROR: ON consensus directory not found: %s\n", consensus_dir_on))
    quit(status = 1)
  }
  if (!dir.exists(consensus_dir_off)) {
    cat(sprintf("ERROR: OFF consensus directory not found: %s\n", consensus_dir_off))
    quit(status = 1)
  }

  available_contrasts <- intersect(
    list.files(consensus_dir_on),
    list.files(consensus_dir_off)
  )

  if (is.null(opts$contrast)) {
    opts$contrast <- available_contrasts[1]
    cat(sprintf("No --contrast specified. Auto-detected: %s\n", opts$contrast))
    cat(sprintf("Available contrasts: %s\n", paste(available_contrasts, collapse = ", ")))
  }

  if (!opts$contrast %in% available_contrasts) {
    cat(sprintf("ERROR: contrast '%s' not found in both ON and OFF directories.\n", opts$contrast))
    cat(sprintf("Available: %s\n", paste(available_contrasts, collapse = ", ")))
    quit(status = 1)
  }

  cat(sprintf("Comparing contrast: %s\n", opts$contrast))
  cat(sprintf("  ON  dir: %s\n", opts$on_dir))
  cat(sprintf("  OFF dir: %s\n", opts$off_dir))
  cat("\n")

  # ── 1. Load consensus results ─────────────────────────────────────────────
  consensus_file_on  <- file.path(consensus_dir_on,  opts$contrast, "consensus_results.tsv")
  consensus_file_off <- file.path(consensus_dir_off, opts$contrast, "consensus_results.tsv")
  summary_file_on    <- file.path(consensus_dir_on,  opts$contrast, "consensus_summary.tsv")
  summary_file_off   <- file.path(consensus_dir_off, opts$contrast, "consensus_summary.tsv")

  cons_on  <- safe_read_tsv(consensus_file_on)
  cons_off <- safe_read_tsv(consensus_file_off)
  sum_on   <- safe_read_tsv(summary_file_on)
  sum_off  <- safe_read_tsv(summary_file_off)

  if (is.null(cons_on) || is.null(cons_off)) {
    cat("ERROR: Could not load consensus results from both directories.\n")
    quit(status = 1)
  }

  cat(sprintf("ON  consensus: %d genes\n", nrow(cons_on)))
  cat(sprintf("OFF consensus: %d genes\n", nrow(cons_off)))

  # ── 2. Consensus summary comparison ───────────────────────────────────────
  tier_counts <- function(df) {
    c(Tier_A       = sum(df$tier == "Tier_A"),
      Tier_B       = sum(df$tier == "Tier_B"),
      Tier_C       = sum(df$tier == "Tier_C"),
      unclassified = sum(df$tier == "unclassified"),
      conflict     = sum(df$consensus_direction == "mixed" | df$consensus_direction == "none"))
  }
  tc_on  <- tier_counts(cons_on)
  tc_off <- tier_counts(cons_off)

  cat("\n── Tier comparison ──\n")
  cat(sprintf("  %-15s %8s %8s %8s\n", "Tier", "ON", "OFF", "Δ"))
  for (tier in names(tc_on)) {
    delta <- tc_off[tier] - tc_on[tier]
    cat(sprintf("  %-15s %8d %8d %+8d\n", tier, tc_on[tier], tc_off[tier], delta))
  }

  # ── 3. Gene-level intersection ────────────────────────────────────────────
  universe <- union(cons_on$gene_id_standard, cons_off$gene_id_standard)
  cat(sprintf("\n── Gene universe ──\n"))
  cat(sprintf("  Union: %d genes\n", length(universe)))

  cons_on_sub  <- cons_on  %>% select(gene_id_standard, consensus_direction, consensus_logFC,
                                       tier, best_rra_fdr, best_cct_fdr, logFC_CV, support_n)
  cons_off_sub <- cons_off %>% select(gene_id_standard, consensus_direction, consensus_logFC,
                                       tier, best_rra_fdr, best_cct_fdr, logFC_CV, support_n)

  colnames(cons_on_sub)  <- paste0(colnames(cons_on_sub),  "_ON")
  colnames(cons_off_sub) <- paste0(colnames(cons_off_sub), "_OFF")

  gene_comp <- full_join(
    cons_on_sub  %>% rename(gene_id_standard = gene_id_standard_ON),
    cons_off_sub %>% rename(gene_id_standard = gene_id_standard_OFF),
    by = "gene_id_standard"
  ) %>%
    mutate(
      logFC_diff      = consensus_logFC_ON - consensus_logFC_OFF,
      logFC_abs_diff  = abs(logFC_diff),
      tier_transition = paste0(tier_ON, "→", tier_OFF),
      in_both         = !is.na(tier_ON) & !is.na(tier_OFF),
      on_only         = !is.na(tier_ON) &  is.na(tier_OFF),
      off_only        =  is.na(tier_ON) & !is.na(tier_OFF),
      # Categorize the transition type
      transition_type = case_when(
        is.na(tier_ON)                              ~ "OFF_only",
        is.na(tier_OFF)                             ~ "ON_only",
        tier_ON == tier_OFF                         ~ "stable",
        tier_ON == "unclassified" & tier_OFF != "unclassified" ~ "promoted_OFF",
        tier_ON != "unclassified" & tier_OFF == "unclassified" ~ "demoted_OFF",
        grepl("Tier", tier_ON) & grepl("Tier", tier_OFF)       ~ "tier_shift",
        TRUE                                       ~ "other"
      )
    )

  n_both    <- sum(gene_comp$in_both, na.rm = TRUE)
  n_on_only <- sum(gene_comp$on_only, na.rm = TRUE)
  n_off_only<- sum(gene_comp$off_only, na.rm = TRUE)

  cat(sprintf("  In both:      %d\n", n_both))
  cat(sprintf("  ON only:      %d\n", n_on_only))
  cat(sprintf("  OFF only:     %d\n", n_off_only))

  # ── 4. Tier transition matrix ─────────────────────────────────────────────
  tier_levels <- c("Tier_A", "Tier_B", "Tier_C", "unclassified")
  tier_mat <- matrix(0, nrow = length(tier_levels) + 1, ncol = length(tier_levels) + 1,
                     dimnames = list(c(tier_levels, "absent_ON"),
                                     c(tier_levels, "absent_OFF")))
  for (i in seq_len(nrow(gene_comp))) {
    rn <- if (is.na(gene_comp$tier_ON[i]))  "absent_ON"  else gene_comp$tier_ON[i]
    cn <- if (is.na(gene_comp$tier_OFF[i])) "absent_OFF" else gene_comp$tier_OFF[i]
    tier_mat[rn, cn] <- tier_mat[rn, cn] + 1
  }

  cat("\n── Tier transition matrix ──\n")
  print(tier_mat)

  tier_mat_df <- as.data.frame(as.table(tier_mat)) %>%
    rename(tier_ON = Var1, tier_OFF = Var2, n_genes = Freq)
  write_tsv(tier_mat_df, file.path(opts$output_dir, "tier_transition_matrix.tsv"))

  # ── 5. Key comparison metrics ─────────────────────────────────────────────
  both_significant <- gene_comp %>%
    filter(in_both, tier_ON != "unclassified", tier_OFF != "unclassified")
  n_both_sig <- nrow(both_significant)

  on_sig_only <- gene_comp %>%
    filter(in_both, tier_ON != "unclassified", tier_OFF == "unclassified")
  n_on_sig_only <- nrow(on_sig_only)

  off_sig_only <- gene_comp %>%
    filter(in_both, tier_ON == "unclassified", tier_OFF != "unclassified")
  n_off_sig_only <- nrow(off_sig_only)

  # Consensus Tier A stability
  on_tier_a  <- gene_comp %>% filter(tier_ON  == "Tier_A") %>% pull(gene_id_standard)
  off_tier_a <- gene_comp %>% filter(tier_OFF == "Tier_A") %>% pull(gene_id_standard)
  tier_a_intersection <- length(intersect(on_tier_a, off_tier_a))
  tier_a_jaccard      <- if (length(union(on_tier_a, off_tier_a)) > 0)
    tier_a_intersection / length(union(on_tier_a, off_tier_a)) else NA_real_
  tier_a_retention    <- if (length(on_tier_a) > 0)
    tier_a_intersection / length(on_tier_a) else NA_real_

  # logFC correlation (genes in both)
  logfc_df <- gene_comp %>%
    filter(in_both, !is.na(consensus_logFC_ON), !is.na(consensus_logFC_OFF))
  logfc_cor <- if (nrow(logfc_df) > 2)
    suppressWarnings(cor(logfc_df$consensus_logFC_ON, logfc_df$consensus_logFC_OFF,
                         use = "complete.obs")) else NA_real_

  cat(sprintf("\n── Key metrics ──\n"))
  cat(sprintf("  Both significant:           %d\n", n_both_sig))
  cat(sprintf("  ON only significant:        %d\n", n_on_sig_only))
  cat(sprintf("  OFF only significant:       %d\n", n_off_sig_only))
  cat(sprintf("  Tier A intersection:        %d\n", tier_a_intersection))
  cat(sprintf("  Tier A Jaccard:             %.4f\n", tier_a_jaccard))
  cat(sprintf("  Tier A retention (ON→OFF):  %.2f%%\n", tier_a_retention * 100))
  cat(sprintf("  logFC Pearson correlation:  %.4f\n", logfc_cor))

  # ── 6. Load decontam efficiency stats ─────────────────────────────────────
  on_stats_dir  <- opts$on_stats_dir  %||% file.path(opts$on_dir,  "03.decontam/stats")
  off_stats_dir <- opts$off_stats_dir %||% file.path(opts$off_dir, "03.decontam/stats")

  # Try fallback: current results/ for ON if backup doesn't have it
  if (!dir.exists(on_stats_dir)) {
    fallback <- opts$stat_dir  # default "results/03.decontam/stats"
    if (dir.exists(fallback)) {
      cat(sprintf("\nNOTE: ON decontam stats not in backup, using %s\n", fallback))
      on_stats_dir <- fallback
    }
  }

  decontam_on  <- collect_decontam_stats(on_stats_dir)
  decontam_off <- collect_decontam_stats(off_stats_dir)

  decontam_efficiency <- NULL
  if (!is.null(decontam_on)) {
    cat("\n── Decontam ON read fate (per-sample) ──\n")
    print(decontam_on %>% select(sample, total_pairs, host_pct, uncertain_pct, removed_pct), n = 20)

    total_summary <- decontam_on %>%
      summarise(
        total_pairs        = sum(total_pairs),
        total_host         = sum(Rescued_Host),
        total_ercc         = sum(Rescued_ERCC),
        total_uncertain    = sum(Flagged_Uncertain),
        total_removed      = sum(Removed_Tech + Removed_NonTarget),
        host_reads_pct     = round(100 * total_host / total_pairs, 2),
        uncertain_pct      = round(100 * total_uncertain / total_pairs, 2),
        removed_pct        = round(100 * total_removed / total_pairs, 2),
        .groups = "drop"
      ) %>% mutate(config = "decontam_ON")

    cat(sprintf("\n  Decontam ON summary:\n"))
    cat(sprintf("    Total pairs:       %d\n", total_summary$total_pairs))
    cat(sprintf("    Host rescued:      %d (%.2f%%)\n", total_summary$total_host, total_summary$host_reads_pct))
    cat(sprintf("    ERCC:              %d\n", total_summary$total_ercc))
    cat(sprintf("    Uncertain:         %d (%.2f%%)\n", total_summary$total_uncertain, total_summary$uncertain_pct))
    cat(sprintf("    Removed:           %d (%.2f%%)\n", total_summary$total_removed, total_summary$removed_pct))

    decontam_efficiency <- decontam_on %>%
      select(sample, total_pairs, Rescued_Host, Rescued_ERCC, Flagged_Uncertain,
             Removed_Tech, Removed_NonTarget, host_pct, uncertain_pct, removed_pct) %>%
      mutate(config = "ON")

    if (!is.null(decontam_off)) {
      # OFF mode: all reads go through (no decontam), but we can still report stats
      # if OFF has its own decontam stats (from a separate analysis), add them
      decontam_off_labeled <- decontam_off %>%
        select(any_of(c("sample", "total_pairs", "Rescued_Host", "Rescued_ERCC",
                        "Flagged_Uncertain", "Removed_Tech", "Removed_NonTarget",
                        "host_pct", "uncertain_pct", "removed_pct"))) %>%
        mutate(config = "OFF")
      decontam_efficiency <- bind_rows(decontam_efficiency, decontam_off_labeled)
    }
  }

  # ── 7. Build comparison summary table ─────────────────────────────────────
  comparison_summary <- tibble(
    metric = c(
      "contrast",
      "gene_universe_union",
      "genes_in_both",
      "genes_ON_only",
      "genes_OFF_only",
      "ON_tier_A",
      "ON_tier_B",
      "ON_tier_C",
      "ON_unclassified",
      "OFF_tier_A",
      "OFF_tier_B",
      "OFF_tier_C",
      "OFF_unclassified",
      "both_significant",
      "ON_only_significant",
      "OFF_only_significant",
      "tier_A_intersection",
      "tier_A_jaccard",
      "tier_A_retention_pct",
      "logFC_pearson_r",
      "decontam_host_retention_pct"
    ),
    value = c(
      opts$contrast,
      length(universe),
      n_both,
      n_on_only,
      n_off_only,
      tc_on["Tier_A"], tc_on["Tier_B"], tc_on["Tier_C"], tc_on["unclassified"],
      tc_off["Tier_A"], tc_off["Tier_B"], tc_off["Tier_C"], tc_off["unclassified"],
      n_both_sig,
      n_on_sig_only,
      n_off_sig_only,
      tier_a_intersection,
      round(tier_a_jaccard, 4),
      round(tier_a_retention * 100, 2),
      round(logfc_cor, 4),
      if (!is.null(decontam_efficiency)) {
        round(100 * sum(decontam_efficiency$Rescued_Host[decontam_efficiency$config == "ON"]) /
                sum(decontam_efficiency$total_pairs[decontam_efficiency$config == "ON"]), 2)
      } else { NA_real_ }
    )
  )

  # ── 8. Write outputs ──────────────────────────────────────────────────────
  write_tsv(comparison_summary,
            file.path(opts$output_dir, "decontam_comparison_summary.tsv"))
  write_tsv(gene_comp,
            file.path(opts$output_dir, "decontam_gene_level_comparison.tsv"))
  if (!is.null(decontam_efficiency)) {
    write_tsv(decontam_efficiency,
              file.path(opts$output_dir, "decontam_efficiency.tsv"))
  }

  # ── 9. Generate figures ───────────────────────────────────────────────────
  fig_dir <- file.path(opts$output_dir, "figures")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  # 9a. logFC scatter: ON vs OFF
  pdf(file.path(fig_dir, "logFC_ON_vs_OFF.pdf"), width = 7, height = 7)
  plot_df <- logfc_df %>%
    mutate(tier_label = case_when(
      tier_ON == "Tier_A" & tier_OFF == "Tier_A" ~ "Tier A (both)",
      tier_ON != "unclassified" & tier_OFF != "unclassified" ~ "Significant (both)",
      tier_ON != "unclassified" & tier_OFF == "unclassified" ~ "ON only",
      tier_ON == "unclassified" & tier_OFF != "unclassified" ~ "OFF only",
      TRUE ~ "Neither"
    ))
  pal <- c("Tier A (both)" = "#2196F3", "Significant (both)" = "#4CAF50",
           "ON only" = "#FF9800", "OFF only" = "#E91E63", "Neither" = "grey80")
  plot(plot_df$consensus_logFC_OFF, plot_df$consensus_logFC_ON,
       col = pal[plot_df$tier_label],
       pch = 16, cex = 0.4,
       xlab = "logFC (decontam OFF)", ylab = "logFC (decontam ON)",
       main = paste0("logFC: decontam ON vs OFF (", opts$contrast, ")"))
  abline(0, 1, lty = 2, col = "grey50")
  abline(h = 0, v = 0, lty = 3, col = "grey70")
  legend("topleft", legend = names(pal), col = pal, pch = 16, cex = 0.7, bty = "n")
  mtext(sprintf("Pearson r = %.3f, n = %d", logfc_cor, nrow(plot_df)),
        side = 3, line = -1.2, cex = 0.8)
  dev.off()
  cat(sprintf("Figure saved: %s\n", file.path(fig_dir, "logFC_ON_vs_OFF.pdf")))

  # 9b. Decontam read fate stacked barplot (if stats available)
  if (!is.null(decontam_on) && nrow(decontam_on) > 0) {
    pdf(file.path(fig_dir, "decontam_read_fate.pdf"), width = max(8, nrow(decontam_on) * 0.8), height = 5)
    fate_mat <- as.matrix(decontam_on[, c("Rescued_Host", "Rescued_ERCC",
                                           "Flagged_Uncertain", "Removed_Tech", "Removed_NonTarget")])
    rownames(fate_mat) <- decontam_on$sample
    fate_prop <- sweep(fate_mat, 1, rowSums(fate_mat), "/")
    barplot(t(fate_prop) * 100, beside = FALSE,
            col = c("#4CAF50", "#2196F3", "#FFC107", "#FF5722", "#9E9E9E"),
            ylab = "Read pairs (%)", xlab = "Sample",
            main = "Decontam read pair fate (ON)",
            las = 2, cex.names = 0.7,
            legend.text = c("Host", "ERCC", "Uncertain", "Tech", "NonTarget"),
            args.legend = list(cex = 0.7, bty = "n"))
    dev.off()
    cat(sprintf("Figure saved: %s\n", file.path(fig_dir, "decontam_read_fate.pdf")))
  }

  # 9c. Tier transition heatmap (text-based matrix plot)
  if (sum(tier_mat) > 0) {
    pdf(file.path(fig_dir, "tier_transition_heatmap.pdf"), width = 6, height = 5.5)
    tier_mat_sub <- tier_mat[c("Tier_A","Tier_B","Tier_C","unclassified","absent_ON"),
                              c("Tier_A","Tier_B","Tier_C","unclassified","absent_OFF")]
    # Normalize by row (ON tier)
    tier_mat_norm <- sweep(tier_mat_sub, 1, rowSums(tier_mat_sub) + 1e-10, "/")
    image(1:ncol(tier_mat_norm), 1:nrow(tier_mat_norm), t(tier_mat_norm[nrow(tier_mat_norm):1, ]),
          col = colorRampPalette(c("white", "#FF5722"))(100),
          axes = FALSE, xlab = "Tier (decontam OFF)", ylab = "Tier (decontam ON)",
          main = paste0("Tier transition: ", opts$contrast))
    axis(1, at = 1:ncol(tier_mat_norm), labels = colnames(tier_mat_norm), las = 2, cex.axis = 0.8)
    axis(2, at = nrow(tier_mat_norm):1, labels = rownames(tier_mat_norm), las = 1, cex.axis = 0.8)
    # Add text labels
    for (ri in 1:nrow(tier_mat_norm)) {
      for (ci in 1:ncol(tier_mat_norm)) {
        text(ci, nrow(tier_mat_norm) - ri + 1,
             labels = tier_mat_sub[ri, ci], cex = 0.7,
             col = if (tier_mat_norm[ri, ci] > 0.5) "white" else "black")
      }
    }
    dev.off()
    cat(sprintf("Figure saved: %s\n", file.path(fig_dir, "tier_transition_heatmap.pdf")))
  }

  # ── 10. Print summary to console ──────────────────────────────────────────
  cat(sprintf("\n── Output files written to %s ──\n", opts$output_dir))
  cat(sprintf("  decontam_comparison_summary.tsv\n"))
  cat(sprintf("  decontam_gene_level_comparison.tsv\n"))
  cat(sprintf("  tier_transition_matrix.tsv\n"))
  if (!is.null(decontam_efficiency))
    cat(sprintf("  decontam_efficiency.tsv\n"))
  cat(sprintf("  figures/logFC_ON_vs_OFF.pdf\n"))
  if (!is.null(decontam_on) && nrow(decontam_on) > 0)
    cat(sprintf("  figures/decontam_read_fate.pdf\n"))
  if (sum(tier_mat) > 0)
    cat(sprintf("  figures/tier_transition_heatmap.pdf\n"))
}

# ══════════════════════════════════════════════════════════════════════════════════
# MODE 4: ALL — run assess, then spike or compare based on result
# ══════════════════════════════════════════════════════════════════════════════════
if (opts$mode == "all") {
  cat("\n=== Decontam Benchmark: FULL PIPELINE ===\n\n")

  # Step 1: Assess
  stat_dir <- opts$stat_dir
  if (!dir.exists(stat_dir)) {
    cat("ERROR: Decontam stats not found. Run pipeline with decontam ON first.\n")
    quit(status = 1)
  }

  classification_files <- list.files(stat_dir, pattern = "_classification_stats.tsv", full.names = TRUE)
  rescue_files         <- list.files(stat_dir, pattern = "_host_rescue_stats.tsv", full.names = TRUE)

  total_host    <- 0
  total_nonhost <- 0
  for (f in classification_files) {
    d <- tryCatch(read_tsv(f, show_col_types = FALSE, comment = "#"), error = function(e) NULL)
    if (is.null(d) || !"unclassified_pairs" %in% names(d)) next
    total_nonhost <- total_nonhost + sum(d$unclassified_pairs, na.rm = TRUE)
  }
  for (f in rescue_files) {
    d <- tryCatch(read_tsv(f, show_col_types = FALSE, comment = "#"), error = function(e) NULL)
    if (is.null(d) || !"host_pairs" %in% names(d)) next
    total_host <- total_host + sum(d$host_pairs, na.rm = TRUE)
  }

  total_pairs <- total_host + total_nonhost
  nonhost_pct <- if (total_pairs > 0) round(100 * total_nonhost / total_pairs, 2) else 0

  cat(sprintf("Non-host reads: %.2f%%\n", nonhost_pct))

  if (nonhost_pct > 5) {
    cat("→ Non-host reads > 5%, using real data comparison.\n\n")
    cat("To complete the comparison, run:\n")
    cat(sprintf("  Rscript workflow/scripts/benchmark_decontam_spikein.R --mode compare \\\n"))
    cat(sprintf("    --on-dir %s --off-dir %s --contrast <contrast> \\\n",
                opts$on_dir, opts$off_dir))
    cat(sprintf("    --output-dir %s\n", opts$output_dir))
  } else {
    if (is.null(opts$microbe_fastq) || !file.exists(opts$microbe_fastq)) {
      cat("→ Non-host reads < 5%, but no microbe FASTQ provided for spike-in.\n")
      cat("  Download Wolbachia reads first:\n")
      cat("    kingfisher get -r SRR1553797 -m ena-ftp -o results/benchmark/\n")
      cat(sprintf("  Then re-run with --microbe-fastq %s/wolbachia.fastq.gz\n", opts$output_dir))
    } else {
      cat("→ Non-host reads < 5%, proceeding with spike-in.\n")
      cat("  (Re-run with --mode spike to execute)\n")
    }
  }
}

cat("\nDone.\n")
