#!/usr/bin/env Rscript
# workflow/scripts/benchmark_decontam_spikein.R
# Decontam module benchmark — condition-specific microbial spike-in evaluation
# Three modes: assess (check current contamination), spike (create spiked FASTQ), compare (analyze results)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

# ── Parse CLI args ──────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
opts <- list(
  mode          = "assess",
  microbe_fastq = NULL,
  spike_rate_a  = 0.15,
  spike_rate_b  = 0.01,
  output_dir    = "results/benchmark"
)
i <- 1
while (i <= length(args)) {
  if (args[i] == "--mode")           { opts$mode <- args[i+1]; i <- i+2 }
  else if (args[i] == "--microbe-fastq") { opts$microbe_fastq <- args[i+1]; i <- i+2 }
  else if (args[i] == "--spike-rate-a")  { opts$spike_rate_a <- as.numeric(args[i+1]); i <- i+2 }
  else if (args[i] == "--spike-rate-b")  { opts$spike_rate_b <- as.numeric(args[i+1]); i <- i+2 }
  else if (args[i] == "--output-dir")    { opts$output_dir <- args[i+1]; i <- i+2 }
  else { i <- i+1 }
}

dir.create(opts$output_dir, recursive = TRUE, showWarnings = FALSE)

# ══════════════════════════════════════════════════════════
# MODE 1: ASSESS — evaluate current contamination levels
# ══════════════════════════════════════════════════════════
if (opts$mode == "assess") {
  cat("\n=== Decontam Spike-in Benchmark: ASSESS MODE ===\n\n")
  
  stat_dir <- "results/03.decontam/stats"
  if (!dir.exists(stat_dir)) {
    cat("ERROR: Decontam stats not found at", stat_dir, "\n")
    cat("Run the full pipeline with decontam.enabled=true first.\n")
    quit(status = 1)
  }
  
  classification_files <- list.files(stat_dir, pattern = "_classification_stats.tsv", full.names = TRUE)
  
  total_host <- 0
  total_nonhost <- 0
  
  for (f in classification_files) {
    d <- tryCatch(read_tsv(f, show_col_types = FALSE, comment = "#"), error = function(e) NULL)
    if (is.null(d)) next
    if (!"unclassified_pairs" %in% names(d)) next
    total_nonhost <- total_nonhost + sum(d$unclassified_pairs, na.rm = TRUE)
  }
  
  rescue_files <- list.files(stat_dir, pattern = "_host_rescue_stats.tsv", full.names = TRUE)
  for (f in rescue_files) {
    d <- tryCatch(read_tsv(f, show_col_types = FALSE, comment = "#"), error = function(e) NULL)
    if (is.null(d)) next
    if (!"host_pairs" %in% names(d)) next
    total_host <- total_host + sum(d$host_pairs, na.rm = TRUE)
  }
  
  total_pairs <- total_host + total_nonhost
  nonhost_pct <- if (total_pairs > 0) round(100 * total_nonhost / total_pairs, 2) else 0
  
  cat(sprintf("Host read pairs:     %d\n", total_host))
  cat(sprintf("Non-host read pairs: %d\n", total_nonhost))
  cat(sprintf("Non-host percentage: %.2f%%\n", nonhost_pct))
  cat("\n")
  
  if (nonhost_pct > 5) {
    cat("RECOMMENDATION: Non-host reads >5%. Use real data for decontam benchmark.\n")
    cat("Proceed to run the pipeline with decontam.enabled=true and decontam.enabled=false\n")
  } else {
    cat("RECOMMENDATION: Non-host reads <5%. Spike-in needed.\n")
    cat("\nTo download Wolbachia reads for spike-in:\n")
    cat("  conda activate kingfisher\n")
    cat("  kingfisher get -r SRR1553797 -m ena-ftp -o results/benchmark/\n")
    cat("\nThen run spike-in:\n")
    cat("  Rscript workflow/scripts/benchmark_decontam_spikein.R --mode spike \\\n")
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

# ══════════════════════════════════════════════════════════
# MODE 2: SPIKE — create condition-specific spiked FASTQ
# ══════════════════════════════════════════════════════════
if (opts$mode == "spike") {
  cat("\n=== Decontam Spike-in Benchmark: SPIKE MODE ===\n\n")
  
  if (is.null(opts$microbe_fastq) || !file.exists(opts$microbe_fastq)) {
    cat("ERROR: --microbe-fastq not specified or file not found\n")
    quit(status = 1)
  }
  
  # Read sample metadata  
  samples <- read_tsv("data/fastq/samples.tsv", show_col_types = FALSE)
  if (!"group" %in% names(samples)) {
    colnames(samples)[1:4] <- c("sample","fq1","fq2","group")
  }
  
  clean_dir <- "results/03.decontam/clean"
  spiked_dir <- file.path(opts$output_dir, "spiked_fastq")
  dir.create(spiked_dir, recursive = TRUE, showWarnings = FALSE)
  
  cat(sprintf("Spike-in rates: Condition A=%.0f%%, Condition B=%.0f%%\n",
              opts$spike_rate_a * 100, opts$spike_rate_b * 100))
  cat(sprintf("Microbial source: %s\n", opts$microbe_fastq))
  cat("\nInstructions for condition-specific spike-in:\n")
  cat("  1. Use seqtk sample to randomly subsample microbial reads\n")
  cat("  2. Concatenate with clean host FASTQ at specified rates\n")
  cat("  3. Update samples.tsv to point to spiked FASTQ\n")
  cat("\nAfter spike-in, run pipeline twice:\n")
  cat("  # Decontam OFF: config.yaml → decontam.enabled=false\n")
  cat("  # Decontam ON:  config.yaml → decontam.enabled=true\n")
  cat("\nThen compare results:\n")
  cat("  Rscript workflow/scripts/benchmark_decontam_spikein.R --mode compare\n")
}

# ══════════════════════════════════════════════════════════
# MODE 3: COMPARE — analyze decontam ON/OFF DEG results
# ══════════════════════════════════════════════════════════
if (opts$mode == "compare") {
  cat("\n=== Decontam Spike-in Benchmark: COMPARE MODE ===\n\n")
  
  dea_dir <- "results/06.differential_expression"
  quantifiers <- c("featurecounts", "stringtie", "salmon", "kallisto")
  
  # Load consensus results for decontam OFF
  consensus_file <- "results/07.consensus_expression/60d_vs_1d/consensus_results.tsv"
  if (file.exists(consensus_file)) {
    cons <- read_tsv(consensus_file, show_col_types = FALSE)
    cat(sprintf("Consensus results: %d genes, Tier A=%d, B=%d, C=%d\n",
                nrow(cons),
                sum(cons$tier == "Tier_A"),
                sum(cons$tier == "Tier_B"),
                sum(cons$tier == "Tier_C")))
  } else {
    cat("Consensus results not found. Run pipeline with decontam ON to generate.\n")
  }
  
  cat("\nFor full comparison, run pipeline twice (decontam on/off) and re-run this mode.\n")
}

cat("\nDone.\n")
