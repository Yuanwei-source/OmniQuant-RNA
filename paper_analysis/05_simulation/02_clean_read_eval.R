#!/usr/bin/env Rscript
library(tximport)
library(DESeq2)
library(RobustRankAggreg)
library(readr)
library(dplyr)

set.seed(42)

BASE <- "experiments/bombyx_enrichment/results/polyester_stress/clean_5rep"
SEEDS <- c(101, 102, 103, 104, 105)
ALPHA <- 0.05
FDR_THRESH <- 0.05
P_CLIP <- 1e-16

tx2gene <- read.table(file.path(BASE, "tx2gene.tsv"), header = FALSE,
                       col.names = c("TXNAME", "GENEID", "GENENAME"),
                       stringsAsFactors = FALSE)
tx2gene <- tx2gene[, c("TXNAME", "GENEID")]

gt_genes_all <- read.delim(file.path(BASE, "ground_truth_genes.tsv"),
                            header = TRUE, stringsAsFactors = FALSE)

condition <- factor(c(rep("control", 5), rep("treatment", 5)))
col_data <- data.frame(condition = condition)
rownames(col_data) <- sprintf("sample_%02d", 1:10)

run_deseq <- function(txi, col_data, name) {
  dds <- DESeqDataSetFromTximport(txi, colData = col_data, design = ~condition)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("condition", "treatment", "control"), alpha = ALPHA)
  res <- as.data.frame(res)
  res$gene_id <- rownames(res)
  res <- res[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat",
                  "pvalue", "padj")]
  rownames(res) <- NULL
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
  res <- res[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat",
                  "pvalue", "padj")]
  rownames(res) <- NULL
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

compute_metrics <- function(predicted, truth, name, decreasing = TRUE) {
  valid <- !is.na(predicted) & !is.na(truth)
  predicted <- predicted[valid]
  truth <- truth[valid]
  
  if (length(predicted) == 0) {
    return(data.frame(method = name, TP = 0, FP = 0, TN = 0, FN = 0,
                      precision = 0, recall = 0, f1 = 0, specificity = 0,
                      n_genes = 0, n_de_true = 0, n_de_called = 0,
                      direction_correct = NA, direction_total = NA,
                      stringsAsFactors = FALSE))
  }
  
  n_de_true <- sum(truth == 1)
  
  if (decreasing) {
    sig <- predicted < ALPHA
  } else {
    sig <- predicted < ALPHA
  }
  
  n_called <- sum(sig)
  
  TP <- sum(sig & truth == 1)
  FP <- sum(sig & truth == 0)
  TN <- sum(!sig & truth == 0)
  FN <- sum(!sig & truth == 1)
  
  precision   <- ifelse(TP + FP > 0, TP / (TP + FP), 0)
  recall      <- ifelse(TP + FN > 0, TP / (TP + FN), 0)
  f1          <- ifelse(precision + recall > 0,
                        2 * precision * recall / (precision + recall), 0)
  specificity <- ifelse(TN + FP > 0, TN / (TN + FP), 0)
  
  data.frame(method = name, TP = TP, FP = FP, TN = TN, FN = FN,
             precision = round(precision, 4), recall = round(recall, 4),
             f1 = round(f1, 4), specificity = round(specificity, 4),
             n_genes = length(truth), n_de_true = n_de_true,
             n_de_called = n_called,
             direction_correct = NA, direction_total = NA,
             stringsAsFactors = FALSE)
}

all_per_method <- list()

for (sim_seed in SEEDS) {
  cat(sprintf("\n========== Seed %d ==========\n", sim_seed))
  
  seed_dir <- file.path(BASE, sprintf("seed_%d", sim_seed))
  
  gt_tx <- read.delim(file.path(seed_dir, "ground_truth.tsv"),
                       header = TRUE, stringsAsFactors = FALSE)
  gt_gene <- aggregate(de_status ~ gene_id, data = gt_tx,
    FUN = function(x) ifelse(any(x == "DE"), "DE", "non-DE"))
  gt_dir  <- aggregate(direction ~ gene_id, data = gt_tx,
    FUN = function(x) unique(x[x != "none"])[1])
  
  gt_genes <- merge(gt_gene, gt_dir, by = "gene_id", all.x = TRUE)
  gt_genes$true_de <- ifelse(gt_genes$de_status == "DE", 1, 0)
  
  cat(sprintf("Gene GT: %d genes, %d DE\n", nrow(gt_genes),
              sum(gt_genes$true_de)))
  
  samples <- sprintf("sample_%02d", 1:10)
  
  salmon_files <- file.path(seed_dir, "salmon_out", samples, "quant.sf")
  names(salmon_files) <- samples
  salmon_txi <- tximport(salmon_files, type = "salmon", tx2gene = tx2gene,
                          countsFromAbundance = "lengthScaledTPM")
  
  kallisto_files <- file.path(seed_dir, "kallisto_out", samples, "abundance.h5")
  if (!file.exists(kallisto_files[1])) {
    kallisto_files <- file.path(seed_dir, "kallisto_out", samples, "abundance.tsv")
  }
  names(kallisto_files) <- samples
  kallisto_txi <- tximport(kallisto_files, type = "kallisto", tx2gene = tx2gene,
                            countsFromAbundance = "lengthScaledTPM")
  
  fc_file <- file.path(seed_dir, "featurecounts_out", "gene_counts.txt")
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
  
  rra_matrix <- cbind(rnk_salmon / N_genes,
                       rnk_kallisto / N_genes,
                       rnk_fc / N_genes)
  rownames(rra_matrix) <- merged$gene_id
  
  rra_result <- aggregateRanks(rmat = rra_matrix, method = "RRA")
  merged$rra_score <- rra_result$Score[match(merged$gene_id, rra_result$Name)]
  
  merged$cct_pvalue <- apply(merged[, c("salmon_pvalue", "kallisto_pvalue",
                                         "featurecounts_pvalue")],
                              1, cauchy_combination)
  merged$cct_pvalue <- as.numeric(merged$cct_pvalue)
  merged$cct_padj <- p.adjust(merged$cct_pvalue, method = "BH")
  
  merged$rra_pseudo_pvalue <- ecdf(merged$rra_score)(merged$rra_score)
  
  fisher_combine <- function(p1, p2) {
    if (is.na(p1) || is.na(p2)) return(NA)
    p1 <- max(p1, 1e-300); p2 <- max(p2, 1e-300)
    chi <- -2 * (log(p1) + log(p2))
    pchisq(chi, df = 4, lower.tail = FALSE)
  }
  merged$fisher_rra_cct <- mapply(fisher_combine,
                                   merged$rra_pseudo_pvalue,
                                   merged$cct_pvalue)
  merged$fisher_rra_cct_padj <- p.adjust(merged$fisher_rra_cct, method = "BH")
  
  merged$tier <- "unclassified"
  merged$rra_sig <- !is.na(merged$rra_score) & merged$rra_score < 0.05
  merged$cct_sig <- !is.na(merged$cct_padj) & merged$cct_padj < 0.05
  merged$fisher_sig <- !is.na(merged$fisher_rra_cct_padj) &
                        merged$fisher_rra_cct_padj < FDR_THRESH
  merged$tier[merged$rra_sig & merged$cct_sig] <- "Tier_A"
  merged$tier[merged$rra_sig & !merged$cct_sig] <- "Tier_B_RRA"
  merged$tier[!merged$rra_sig & merged$cct_sig] <- "Tier_B_CCT"
  
  logfc_cols <- c("salmon_log2FC", "kallisto_log2FC", "featurecounts_log2FC")
  merged$n_pos <- rowSums(!is.na(merged[, logfc_cols]) &
                           merged[, logfc_cols] > 0)
  merged$n_neg <- rowSums(!is.na(merged[, logfc_cols]) &
                           merged[, logfc_cols] < 0)
  merged$consensus_log2FC <- rowMeans(merged[, logfc_cols], na.rm = TRUE)
  merged$sign_consistent <- merged$n_pos == 3 | merged$n_neg == 3
  
  eval_data <- merge(merged, gt_genes, by = "gene_id", all.x = TRUE)
  
  methods <- list(
    list(col = "salmon_pvalue",       name = "Salmon",       dec = TRUE),
    list(col = "kallisto_pvalue",     name = "Kallisto",     dec = TRUE),
    list(col = "featurecounts_pvalue", name = "featureCounts", dec = TRUE),
    list(col = "cct_pvalue",          name = "CCT",          dec = TRUE),
    list(col = "rra_score",           name = "RRA",          dec = FALSE),
    list(col = "fisher_rra_cct",      name = "Fisher(RRA+CCT)", dec = TRUE)
  )
  
  seed_metrics <- list()
  for (m in methods) {
    mt <- compute_metrics(eval_data[[m$col]], eval_data$true_de, m$name,
                          decreasing = m$dec)
    seed_metrics[[m$name]] <- mt
    cat(sprintf("  %-18s: P=%.4f  R=%.4f  F1=%.4f  called=%d/%d\n",
      m$name, mt$precision, mt$recall, mt$f1, mt$n_de_called, mt$n_de_true))
  }
  
  ta_sig <- eval_data$tier == "Tier_A"
  ta_tp <- sum(ta_sig & eval_data$true_de == 1, na.rm = TRUE)
  ta_fp <- sum(ta_sig & eval_data$true_de == 0, na.rm = TRUE)
  ta_tn <- sum(!ta_sig & eval_data$true_de == 0, na.rm = TRUE)
  ta_fn <- sum(!ta_sig & eval_data$true_de == 1, na.rm = TRUE)
  ta_prec <- ifelse(ta_tp + ta_fp > 0, ta_tp / (ta_tp + ta_fp), 0)
  ta_rec  <- ifelse(ta_tp + ta_fn > 0, ta_tp / (ta_tp + ta_fn), 0)
  ta_f1   <- ifelse(ta_prec + ta_rec > 0,
                     2 * ta_prec * ta_rec / (ta_prec + ta_rec), 0)
  
  tier_a_metrics <- data.frame(
    method = "RRA+CCT(TierA)", TP = ta_tp, FP = ta_fp, TN = ta_tn, FN = ta_fn,
    precision = round(ta_prec, 4), recall = round(ta_rec, 4),
    f1 = round(ta_f1, 4), specificity = round(ta_tn / (ta_tn + ta_fp), 4),
    n_genes = nrow(eval_data), n_de_true = sum(eval_data$true_de == 1),
    n_de_called = sum(ta_sig, na.rm = TRUE),
    direction_correct = NA, direction_total = NA,
    stringsAsFactors = FALSE)
  seed_metrics[["RRA+CCT(TierA)"]] <- tier_a_metrics
  cat(sprintf("  %-18s: P=%.4f  R=%.4f  F1=%.4f  called=%d/%d (Tier A)\n",
    "RRA+CCT(TierA)", ta_prec, ta_rec, ta_f1,
    sum(ta_sig, na.rm = TRUE), sum(eval_data$true_de == 1, na.rm = TRUE)))
  
  mdf <- do.call(rbind, seed_metrics)
  mdf$seed <- sim_seed
  all_per_method[[as.character(sim_seed)]] <- mdf
  
  out_dir <- file.path(seed_dir, "dea_results")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  write.table(mdf, file.path(out_dir, "metrics.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(eval_data, file.path(out_dir, "full_results.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

cat("\n========== AGGREGATE RESULTS ==========\n")
agg <- do.call(rbind, all_per_method)

for (m in unique(agg$method)) {
  mm <- agg[agg$method == m, ]
  cat(sprintf("\n--- %s (over %d seeds) ---\n", m, nrow(mm)))
  cat(sprintf("  F1:      %.4f ± %.4f\n", mean(mm$f1), sd(mm$f1)))
  cat(sprintf("  Prec:    %.4f ± %.4f\n", mean(mm$precision), sd(mm$precision)))
  cat(sprintf("  Recall:  %.4f ± %.4f\n", mean(mm$recall), sd(mm$recall)))
  cat(sprintf("  Called:  %.0f / %.0f true DE\n", mean(mm$n_de_called),
              mean(mm$n_de_true)))
}

summary_dir <- file.path(BASE, "summary")
dir.create(summary_dir, showWarnings = FALSE, recursive = TRUE)

write.table(agg, file.path(summary_dir, "all_metrics.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

summary_table <- data.frame(
  Method = unique(agg$method),
  stringsAsFactors = FALSE
)
summary_table$Precision_mean <- sapply(summary_table$Method, function(m) {
  round(mean(agg$precision[agg$method == m]), 4)
})
summary_table$Precision_sd <- sapply(summary_table$Method, function(m) {
  round(sd(agg$precision[agg$method == m]), 4)
})
summary_table$Recall_mean <- sapply(summary_table$Method, function(m) {
  round(mean(agg$recall[agg$method == m]), 4)
})
summary_table$Recall_sd <- sapply(summary_table$Method, function(m) {
  round(sd(agg$recall[agg$method == m]), 4)
})
summary_table$F1_mean <- sapply(summary_table$Method, function(m) {
  round(mean(agg$f1[agg$method == m]), 4)
})
summary_table$F1_sd <- sapply(summary_table$Method, function(m) {
  round(sd(agg$f1[agg$method == m]), 4)
})
summary_table$Called_mean <- sapply(summary_table$Method, function(m) {
  round(mean(agg$n_de_called[agg$method == m]), 1)
})
summary_table$Called_sd <- sapply(summary_table$Method, function(m) {
  round(sd(agg$n_de_called[agg$method == m]), 1)
})
summary_table$Direction <- c(
  "bidirectional", "bidirectional", "bidirectional",
  "bidirectional", "bidirectional", "bidirectional", "bidirectional"
)

write.table(summary_table, file.path(summary_dir, "summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n========== FINAL SUMMARY TABLE ==========\n")
cat(sprintf("%-20s %14s %14s %14s %14s\n",
            "Method", "Precision", "Recall", "F1", "Called"))
cat(strrep("-", 80), "\n")
for (i in seq_len(nrow(summary_table))) {
  cat(sprintf("%-20s %5.4f(%.4f) %5.4f(%.4f) %5.4f(%.4f) %5.0f(%.1f)\n",
    summary_table$Method[i],
    summary_table$Precision_mean[i], summary_table$Precision_sd[i],
    summary_table$Recall_mean[i], summary_table$Recall_sd[i],
    summary_table$F1_mean[i], summary_table$F1_sd[i],
    summary_table$Called_mean[i], summary_table$Called_sd[i]))
}

cat(sprintf("\nResults saved to: %s\n", summary_dir))
cat("=== B-lite benchmark evaluation complete ===\n")
