#!/usr/bin/env Rscript
# =============================================================================
# Stress Test 5-repeat: DEA + consensus evaluation
# Per seed: DESeq2 + RRA + CCT + Tier A metrics
# Aggregate: mean ± sd across 5 seeds
# BUG FIX: gene_id as FIRST column in DESeq2 output
# =============================================================================
suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(RobustRankAggreg)
  library(readr)
  library(dplyr)
})
set.seed(42)

BASE <- "experiments/bombyx_enrichment/results/polyester_stress/stress_5rep"
SEEDS <- c(201, 202, 203, 204, 205)
ALPHA <- 0.05
FDR_THRESH <- 0.05
N_SAMPLES <- 10

tx2gene <- read.table(file.path(BASE, "ground_truth_genes.tsv"),
                       header = TRUE, stringsAsFactors = FALSE)

tx2gene_map <- read.table(
  "experiments/bombyx_enrichment/results/polyester_stress/stress_v1/tx2gene.tsv",
  header = FALSE, col.names = c("tx_id", "tx_name", "gene_id"),
  stringsAsFactors = FALSE)
tx2gene_sub <- tx2gene_map[, c("tx_id", "gene_id")]

condition <- factor(c(rep("control", 5), rep("treatment", 5)))
col_data <- data.frame(condition = condition)
rownames(col_data) <- sprintf("sample_%02d", 1:N_SAMPLES)

run_deseq <- function(txi, col_data, name) {
  dds <- DESeqDataSetFromTximport(txi, colData = col_data, design = ~condition)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("condition", "treatment", "control"), alpha = ALPHA)
  res <- as.data.frame(res)
  res$gene_id <- rownames(res)
  rownames(res) <- NULL
  res <- res[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
  colnames(res) <- c("gene_id", "baseMean",
                     paste0(name, "_log2FC"), paste0(name, "_lfcSE"),
                     paste0(name, "_stat"), paste0(name, "_pvalue"),
                     paste0(name, "_padj"))
  res
}

run_deseq_counts <- function(count_matrix, col_data, name) {
  dds <- DESeqDataSetFromMatrix(countData = round(count_matrix),
                                colData = col_data, design = ~condition)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("condition", "treatment", "control"), alpha = ALPHA)
  res <- as.data.frame(res)
  res$gene_id <- rownames(res)
  rownames(res) <- NULL
  res <- res[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
  colnames(res) <- c("gene_id", "baseMean",
                     paste0(name, "_log2FC"), paste0(name, "_lfcSE"),
                     paste0(name, "_stat"), paste0(name, "_pvalue"),
                     paste0(name, "_padj"))
  res
}

cauchy_combination <- function(pvals) {
  pvals <- as.numeric(pvals)
  valid <- !is.na(pvals) & is.finite(pvals)
  if (sum(valid) == 0) return(NA)
  pvals <- pvals[valid]
  pvals[pvals < 1e-300] <- 1e-300
  pvals[pvals > (1 - 1e-15)] <- 1 - 1e-15
  weights <- rep(1 / length(pvals), length(pvals))
  cauchy_vals <- tan(pi * (0.5 - pvals))
  stat <- sum(weights * cauchy_vals)
  pval <- 1 - pcauchy(stat)
  if (pval > 0.5) pval <- 1 - pval
  pval <- 2 * min(pval, 1 - pval)
  max(pval, 1e-300)
}

fisher_combine <- function(p1, p2) {
  if (is.na(p1) || is.na(p2)) return(NA)
  p1 <- max(p1, 1e-300); p2 <- max(p2, 1e-300)
  chi <- -2 * (log(p1) + log(p2))
  pchisq(chi, df = 4, lower.tail = FALSE)
}

compute_metrics <- function(pval_vec, true_de_vec, name) {
  valid <- !is.na(pval_vec) & !is.na(true_de_vec)
  pval_vec <- pval_vec[valid]
  true_de_vec <- true_de_vec[valid]
  
  sig <- pval_vec < FDR_THRESH
  n_called <- sum(sig)
  
  TP <- sum(sig & true_de_vec == 1)
  FP <- sum(sig & true_de_vec == 0)
  TN <- sum(!sig & true_de_vec == 0)
  FN <- sum(!sig & true_de_vec == 1)
  
  precision <- ifelse(TP + FP > 0, TP / (TP + FP), 0)
  recall    <- ifelse(TP + FN > 0, TP / (TP + FN), 0)
  f1        <- ifelse(precision + recall > 0, 2 * precision * recall / (precision + recall), 0)
  
  list(name = name, TP = TP, FP = FP, TN = TN, FN = FN,
       precision = precision, recall = recall, f1 = f1, called = n_called)
}

compute_iso_switch_fp <- function(pval_vec, gene_nums) {
  iso_switch <- gene_nums >= 101 & gene_nums <= 200
  valid <- !is.na(pval_vec) & iso_switch
  sig <- pval_vec[valid] < FDR_THRESH
  fp_count <- sum(sig, na.rm = TRUE)
  total <- sum(valid, na.rm = TRUE)
  list(count = fp_count, total = total, rate = ifelse(total > 0, fp_count / total, 0))
}

all_per_method <- list()

for (sim_seed in SEEDS) {
  cat(sprintf("\n========== Seed %d ==========\n", sim_seed))
  
  seed_dir <- file.path(BASE, sprintf("seed_%d", sim_seed))
  
  gt_genes <- read.table(file.path(seed_dir, "ground_truth_genes.tsv"),
                          header = TRUE, stringsAsFactors = FALSE)
  
  cat(sprintf("Gene GT: %d genes, %d true DE\n", nrow(gt_genes), sum(gt_genes$true_de)))
  
  samples <- sprintf("sample_%02d", 1:N_SAMPLES)
  
  salmon_files <- file.path(seed_dir, "salmon_out", samples, "quant.sf")
  names(salmon_files) <- samples
  if (!all(file.exists(salmon_files))) {
    cat("WARNING: Missing salmon files, skipping seed\n")
    next
  }
  salmon_txi <- tximport(salmon_files, type = "salmon", tx2gene = tx2gene_sub,
                         ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
  
  kallisto_files <- file.path(seed_dir, "kallisto_out", samples, "abundance.tsv")
  names(kallisto_files) <- samples
  if (!all(file.exists(kallisto_files))) {
    cat("WARNING: Missing kallisto files, skipping seed\n")
    next
  }
  kallisto_txi <- tximport(kallisto_files, type = "kallisto", tx2gene = tx2gene_sub,
                            ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
  
  fc_file <- file.path(seed_dir, "featurecounts_out", "gene_counts.txt")
  if (!file.exists(fc_file)) {
    cat("WARNING: Missing featureCounts file, skipping seed\n")
    next
  }
  fc_raw <- read.delim(fc_file, header = TRUE, comment.char = "#",
                        stringsAsFactors = FALSE, check.names = FALSE)
  rownames(fc_raw) <- fc_raw[, 1]
  count_cols <- grep("sorted.bam$", colnames(fc_raw), value = TRUE)
  fc_counts <- as.matrix(fc_raw[, count_cols, drop = FALSE])
  colnames(fc_counts) <- samples
  rownames(fc_counts) <- fc_raw[, 1]
  
  salmon_res    <- run_deseq(salmon_txi, col_data, "salmon")
  kallisto_res  <- run_deseq(kallisto_txi, col_data, "kallisto")
  fc_res        <- run_deseq_counts(fc_counts, col_data, "featurecounts")
  
  cat(sprintf("DESeq2 sig(p<0.05): salmon=%d, kallisto=%d, fc=%d\n",
    sum(salmon_res$salmon_pvalue < 0.05, na.rm = TRUE),
    sum(kallisto_res$kallisto_pvalue < 0.05, na.rm = TRUE),
    sum(fc_res$featurecounts_pvalue < 0.05, na.rm = TRUE)))
  
  merged <- merge(salmon_res, kallisto_res, by = "gene_id", all = TRUE)
  merged <- merge(merged, fc_res, by = "gene_id", all = TRUE)
  
  N_genes <- nrow(merged)
  rnk_salmon <- rank(merged$salmon_pvalue, ties.method = "min", na.last = "keep")
  rnk_kallisto <- rank(merged$kallisto_pvalue, ties.method = "min", na.last = "keep")
  rnk_fc <- rank(merged$featurecounts_pvalue, ties.method = "min", na.last = "keep")
  
  rra_matrix <- cbind(rnk_salmon / N_genes, rnk_kallisto / N_genes, rnk_fc / N_genes)
  rownames(rra_matrix) <- merged$gene_id
  
  rra_result <- aggregateRanks(rmat = rra_matrix, method = "RRA")
  merged$rra_score <- rra_result$Score[match(merged$gene_id, rra_result$Name)]
  
  merged$cct_pvalue <- apply(merged[, c("salmon_pvalue", "kallisto_pvalue", "featurecounts_pvalue")],
                              1, cauchy_combination)
  merged$cct_pvalue <- as.numeric(merged$cct_pvalue)
  merged$cct_padj <- p.adjust(merged$cct_pvalue, method = "BH")
  
  merged$rra_pseudo_pvalue <- ecdf(merged$rra_score)(merged$rra_score)
  merged$fisher_rra_cct <- mapply(fisher_combine, merged$rra_pseudo_pvalue, merged$cct_pvalue)
  merged$fisher_rra_cct_padj <- p.adjust(merged$fisher_rra_cct, method = "BH")
  
  merged$rra_sig <- !is.na(merged$rra_score) & merged$rra_score < 0.05
  merged$cct_sig <- !is.na(merged$cct_padj) & merged$cct_padj < FDR_THRESH
  merged$fisher_sig <- !is.na(merged$fisher_rra_cct_padj) & merged$fisher_rra_cct_padj < FDR_THRESH
  
  merged$tier <- "unclassified"
  merged$tier[merged$rra_sig & merged$cct_sig] <- "Tier_A"
  merged$tier[merged$rra_sig & !merged$cct_sig] <- "Tier_B_RRA"
  merged$tier[!merged$rra_sig & merged$cct_sig] <- "Tier_B_CCT"
  
  eval_data <- merge(merged, gt_genes, by = "gene_id", all.x = TRUE)
  
  gene_num_map <- setNames(gt_genes$gene_num, gt_genes$gene_id)
  eval_data$gene_num <- gene_num_map[eval_data$gene_id]
  
  methods <- list(
    list(col = "salmon_pvalue",       name = "Salmon"),
    list(col = "kallisto_pvalue",     name = "Kallisto"),
    list(col = "featurecounts_pvalue", name = "featureCounts"),
    list(col = "fisher_rra_cct",      name = "Fisher(RRA+CCT)")
  )
  
  seed_metrics <- list()
  for (m in methods) {
    res <- compute_metrics(eval_data[[m$col]], eval_data$true_de, m$name)
    iso_fp <- compute_iso_switch_fp(eval_data[[m$col]], eval_data$gene_num)
    res$iso_fp_count <- iso_fp$count
    res$iso_fp_total <- iso_fp$total
    res$iso_fp_rate <- iso_fp$rate
    seed_metrics[[m$name]] <- res
    cat(sprintf("  %-18s: P=%.4f R=%.4f F1=%.4f IsoFP=%.0f%%(%d/%d) called=%d\n",
      m$name, res$precision, res$recall, res$f1,
      res$iso_fp_rate * 100, res$iso_fp_count, res$iso_fp_total,
      res$called))
  }
  
  ta_sig <- eval_data$tier == "Tier_A"
  ta_tp <- sum(ta_sig & eval_data$true_de == 1, na.rm = TRUE)
  ta_fp <- sum(ta_sig & eval_data$true_de == 0, na.rm = TRUE)
  ta_tn <- sum(!ta_sig & eval_data$true_de == 0, na.rm = TRUE)
  ta_fn <- sum(!ta_sig & eval_data$true_de == 1, na.rm = TRUE)
  ta_prec <- ifelse(ta_tp + ta_fp > 0, ta_tp / (ta_tp + ta_fp), 0)
  ta_rec  <- ifelse(ta_tp + ta_fn > 0, ta_tp / (ta_tp + ta_fn), 0)
  ta_f1   <- ifelse(ta_prec + ta_rec > 0, 2 * ta_prec * ta_rec / (ta_prec + ta_rec), 0)
  
  ta_iso_switch <- eval_data$gene_num >= 101 & eval_data$gene_num <= 200
  ta_iso_fp <- sum(ta_sig & ta_iso_switch, na.rm = TRUE)
  ta_iso_total <- sum(ta_iso_switch, na.rm = TRUE)
  
  tier_a_res <- list(
    name = "RRA+CCT(TierA)", TP = ta_tp, FP = ta_fp, TN = ta_tn, FN = ta_fn,
    precision = ta_prec, recall = ta_rec, f1 = ta_f1,
    called = sum(ta_sig, na.rm = TRUE),
    iso_fp_count = ta_iso_fp, iso_fp_total = ta_iso_total,
    iso_fp_rate = ifelse(ta_iso_total > 0, ta_iso_fp / ta_iso_total, 0)
  )
  seed_metrics[["RRA+CCT(TierA)"]] <- tier_a_res
  
  cat(sprintf("  %-18s: P=%.4f R=%.4f F1=%.4f IsoFP=%.0f%%(%d/%d) called=%d (Tier A)\n",
    "RRA+CCT(TierA)", ta_prec, ta_rec, ta_f1,
    tier_a_res$iso_fp_rate * 100, ta_iso_fp, ta_iso_total,
    tier_a_res$called))
  
  mdf <- do.call(rbind, lapply(seed_metrics, function(r) {
    data.frame(method = r$name,
               precision = r$precision, recall = r$recall, f1 = r$f1,
               iso_fp_rate = r$iso_fp_rate, iso_fp_count = r$iso_fp_count,
               iso_fp_total = r$iso_fp_total, called = r$called,
               TP = r$TP, FP = r$FP, stringsAsFactors = FALSE)
  }))
  mdf$seed <- sim_seed
  all_per_method[[as.character(sim_seed)]] <- mdf
  
  out_dir <- file.path(seed_dir, "dea_results")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  write.table(mdf, file.path(out_dir, "metrics.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(eval_data, file.path(out_dir, "full_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}

cat("\n========== AGGREGATE RESULTS ==========\n")
agg <- do.call(rbind, all_per_method)

display_order <- c("Salmon", "Kallisto", "featureCounts", "Fisher(RRA+CCT)", "RRA+CCT(TierA)")

cat("\n=== Per-method summary (mean ± sd over 5 seeds) ===\n\n")
cat(sprintf("%-20s %14s %14s %14s %14s\n",
            "Method", "Precision", "Recall", "F1", "IsoFP_pct"))
cat(strrep("-", 80), "\n")

for (nm in display_order) {
  mm <- agg[agg$method == nm, ]
  if (nrow(mm) == 0) next
  p_mean <- mean(mm$precision); p_sd <- sd(mm$precision)
  r_mean <- mean(mm$recall);    r_sd <- sd(mm$recall)
  f_mean <- mean(mm$f1);        f_sd <- sd(mm$f1)
  ifp_mean <- mean(mm$iso_fp_rate * 100); ifp_sd <- sd(mm$iso_fp_rate * 100)
  
  cat(sprintf("%-20s %5.3f(%5.3f) %5.3f(%5.3f) %5.3f(%5.3f) %5.1f(%5.1f)\n",
    nm, p_mean, p_sd, r_mean, r_sd, f_mean, f_sd, ifp_mean, ifp_sd))
}

cat("\n=== Called genes (mean ± sd) ===\n")
for (nm in display_order) {
  mm <- agg[agg$method == nm, ]
  if (nrow(mm) == 0) next
  cat(sprintf("  %-20s: %.0f ± %.0f\n", nm, mean(mm$called), sd(mm$called)))
}

cat("\n=== Iso-switch FP (mean ± sd) ===\n")
for (nm in display_order) {
  mm <- agg[agg$method == nm, ]
  if (nrow(mm) == 0) next
  cat(sprintf("  %-20s: %.1f%% ± %.1f%% (%d/%d avg)\n", nm,
    mean(mm$iso_fp_rate * 100), sd(mm$iso_fp_rate * 100),
    as.integer(mean(mm$iso_fp_count)), as.integer(mean(mm$iso_fp_total))))
}

cat("\n=== Tier A vs best single quantifier ===\n")
ta_mm <- agg[agg$method == "RRA+CCT(TierA)", ]
# Find best single quantifier by F1
single_methods <- c("Salmon", "Kallisto", "featureCounts")
best_single <- single_methods[1]
best_single_f1 <- 0
for (nm in single_methods) {
  sm <- agg[agg$method == nm, ]
  if (nrow(sm) > 0 && mean(sm$f1) > best_single_f1) {
    best_single_f1 <- mean(sm$f1)
    best_single <- nm
  }
}

if (nrow(ta_mm) > 0) {
  bs_mm <- agg[agg$method == best_single, ]
  
  delta_prec <- mean(ta_mm$precision) - mean(bs_mm$precision)
  delta_f1   <- mean(ta_mm$f1) - mean(bs_mm$f1)
  delta_ifp  <- mean(bs_mm$iso_fp_rate * 100) - mean(ta_mm$iso_fp_rate * 100)
  
  ta_prec_dir <- all(ta_mm$precision > bs_mm$precision)
  ta_ifp_dir  <- all(ta_mm$iso_fp_rate < bs_mm$iso_fp_rate)
  
  cat(sprintf("  Best single quantifier: %s (F1=%.3f)\n", best_single, best_single_f1))
  cat(sprintf("  ΔPrecision (TierA - %s): %+.3f pp (%d/5 seeds direction consistent: %s)\n",
    best_single, delta_prec * 100, nrow(ta_mm), ifelse(ta_prec_dir, "YES", "NO")))
  cat(sprintf("  ΔF1 (TierA - %s): %+.3f\n", best_single, delta_f1))
  cat(sprintf("  ΔIsoFP (%s - TierA): %+.1f pp (%d/5 seeds direction consistent: %s)\n",
    best_single, delta_ifp, nrow(ta_mm), ifelse(ta_ifp_dir, "YES", "NO")))
  
  cat(sprintf("\n  Tier A mean IsoFP: %.1f%% ± %.1f%%\n",
    mean(ta_mm$iso_fp_rate * 100), sd(ta_mm$iso_fp_rate * 100)))
  cat(sprintf("  %s mean IsoFP: %.1f%% ± %.1f%%\n",
    best_single, mean(bs_mm$iso_fp_rate * 100), sd(bs_mm$iso_fp_rate * 100)))
}

cat("\n=== Verdict ===\n")
if (nrow(ta_mm) > 0 && nrow(bs_mm) > 0) {
  if (mean(ta_mm$iso_fp_rate) < mean(bs_mm$iso_fp_rate)) {
    cat(sprintf("PASS - Tier A iso-switch FP%% (%.0f%%) is lower than %s (%.0f%%)\n",
      mean(ta_mm$iso_fp_rate) * 100, best_single, mean(bs_mm$iso_fp_rate) * 100))
  } else {
    cat(sprintf("WARNING - Tier A iso-switch FP%% (%.0f%%) NOT lower than %s (%.0f%%)\n",
      mean(ta_mm$iso_fp_rate) * 100, best_single, mean(bs_mm$iso_fp_rate) * 100))
  }
}

summary_dir <- file.path(BASE, "summary")
dir.create(summary_dir, showWarnings = FALSE, recursive = TRUE)
write.table(agg, file.path(summary_dir, "all_metrics.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

summary_table <- data.frame(
  Method = display_order,
  stringsAsFactors = FALSE
)
for (col_name in c("precision", "recall", "f1", "iso_fp_rate")) {
  summary_table[[paste0(col_name, "_mean")]] <- sapply(summary_table$Method, function(m) {
    mm <- agg[agg$method == m, ]
    if (nrow(mm) == 0) return(NA)
    round(mean(mm[[col_name]]), 4)
  })
  summary_table[[paste0(col_name, "_sd")]] <- sapply(summary_table$Method, function(m) {
    mm <- agg[agg$method == m, ]
    if (nrow(mm) == 0) return(NA)
    round(sd(mm[[col_name]]), 4)
  })
}
write.table(summary_table, file.path(summary_dir, "summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("\nResults saved to: %s\n", summary_dir))
cat("=== Stress test 5-repeat evaluation complete ===\n")
