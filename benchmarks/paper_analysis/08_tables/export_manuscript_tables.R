#!/usr/bin/env Rscript
# Export all manuscript tables from final summary files
suppressPackageStartupMessages({ library(dplyr); library(readr) })

OUT <- "paper_outputs/tables_main"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

# Table: annotation degradation
cat("annotation degradation: see benchmark report\n")

# Table: decontam ON/OFF
cat("decontam ON/OFF: see benchmark report\n")

# Table: clean benchmark
clean <- tryCatch(read_tsv("experiments/bombyx_enrichment/results/polyester_stress/clean_5rep/summary/summary.tsv", show_col_types=FALSE), error=function(e)NULL)
if(!is.null(clean)) write_tsv(clean, file.path(OUT, "table_clean_benchmark.tsv"))

# Table: stress benchmark
stress <- tryCatch(read_tsv("experiments/bombyx_enrichment/results/polyester_stress/stress_5rep/summary/summary.tsv", show_col_types=FALSE), error=function(e)NULL)
if(!is.null(stress)) write_tsv(stress, file.path(OUT, "table_stress_benchmark.tsv"))

# Table: SEQC validation
seqc <- tryCatch(
  read_tsv("experiments/bombyx_enrichment/results/simulated_benchmark/seqc_validation.tsv", show_col_types=FALSE),
  error = function(e) {
    cat(sprintf("  WARNING: SEQC validation file not found, skipping\n"))
    NULL
  }
)
if(!is.null(seqc)) write_tsv(seqc, file.path(OUT, "table_seqc_validation.tsv"))

cat("Tables exported.\n")
