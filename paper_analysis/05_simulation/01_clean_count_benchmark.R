#!/usr/bin/env Rscript
###############################################################################
# Phase 1A/B: Simulated count data benchmark for OmniQuant-RNA
#
# Compares DESeq2, edgeR, limma-voom individually vs RRA/CCT/Fisher/Stouffer
# consensus methods on simulated count data with known ground truth.
#
# - Null simulations (0% DEG): Type I error calibration (500 repeats)
# - Non-null simulations (5/10/20% DEG): PR curves, F1, direction concordance
###############################################################################

suppressPackageStartupMessages({
  library(DESeq2); library(edgeR); library(limma); library(RobustRankAggreg)
  library(dplyr); library(readr); library(tidyr); library(tibble)
})

set.seed(42)

# ── Configuration ────────────────────────────────────────────────────────────
OUT_DIR <- "experiments/bombyx_enrichment/results/simulated_benchmark"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

N_GENES <- 10000L
N_REP   <- 3L
P_CLIP  <- 1e-16
ALPHA   <- 0.05

# ── Count simulation from negative binomial ──────────────────────────────────
simulate_counts <- function(n_genes = N_GENES, n_rep = N_REP, deg_pct = 0,
                            log2fc = 1.0, base_mean = 1000, disp = 0.1) {
  n_de <- round(n_genes * deg_pct)
  n_samples <- 2 * n_rep
  mu <- matrix(rgamma(n_genes, shape = 20, rate = 20/base_mean), n_genes, 1)
  counts <- matrix(NA_integer_, n_genes, n_samples)
  for (i in seq_len(n_genes)) {
    for (j in seq_len(n_samples)) {
      counts[i, j] <- rnbinom(1, mu = mu[i], size = 1/disp)
    }
  }
  truth <- rep(0L, n_genes)
  if (n_de > 0) {
    de_idx   <- sample(n_genes, n_de, replace = FALSE)
    half     <- floor(n_de / 2)
    up_idx   <- de_idx[1:half]
    down_idx <- de_idx[(half + 1):n_de]
    for (i in up_idx) {
      counts[i, (n_rep + 1):n_samples] <-
        rnbinom(n_rep, mu = mu[i] * 2^log2fc, size = 1/disp)
    }
    for (i in down_idx) {
      counts[i, (n_rep + 1):n_samples] <-
        rnbinom(n_rep, mu = mu[i] * 2^(-log2fc), size = 1/disp)
    }
    truth[up_idx]   <- 1L
    truth[down_idx] <- -1L
  }
  colnames(counts) <- c(paste0("Ctrl_", 1:n_rep), paste0("Trt_", 1:n_rep))
  rownames(counts) <- paste0("gene_", 1:n_genes)
  list(counts = counts, truth = truth, de_idx = which(truth != 0))
}

# ── Single method DEA ───────────────────────────────────────────────────────
run_deseq2 <- function(counts, n_rep) {
  coldata <- data.frame(condition = factor(rep(c("Ctrl", "Trt"), each = n_rep)))
  dds <- DESeqDataSetFromMatrix(round(counts), coldata, ~condition)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("condition", "Trt", "Ctrl"))
  data.frame(gene = rownames(res), logFC = res$log2FoldChange,
             p = res$pvalue, padj = res$padj, stringsAsFactors = FALSE)
}

run_edger <- function(counts, n_rep) {
  y <- DGEList(round(counts), group = rep(c("Ctrl","Trt"), each = n_rep))
  y <- calcNormFactors(y)
  design <- model.matrix(~group, data = data.frame(group = y$samples$group))
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = 2)
  data.frame(gene = rownames(qlf$table), logFC = qlf$table$logFC,
             p = qlf$table$PValue, padj = p.adjust(qlf$table$PValue, "BH"),
             stringsAsFactors = FALSE)
}

run_limma <- function(counts, n_rep) {
  dge <- DGEList(round(counts), group = rep(c("Ctrl","Trt"), each = n_rep))
  dge <- calcNormFactors(dge)
  design <- model.matrix(~rep(c("Ctrl","Trt"), each = n_rep))
  v <- voom(dge, design)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef = 2, number = nrow(counts), sort.by = "none")
  data.frame(gene = rownames(tt), logFC = tt$logFC, p = tt$P.Value,
             padj = tt$adj.P.Val, stringsAsFactors = FALSE)
}

# ── P-value helpers ─────────────────────────────────────────────────────────
p_one <- function(p_val, logfc, direction = "up") {
  if (direction == "up") ifelse(logfc > 0, p_val/2, 1 - p_val/2)
  else ifelse(logfc < 0, p_val/2, 1 - p_val/2)
}

fisher_combine <- function(pvals) {
  statistic <- -2 * sum(log(pmax(pvals, .Machine$double.xmin)))
  pchisq(statistic, 2 * length(pvals), lower.tail = FALSE)
}

stouffer_combine <- function(pvals) {
  z <- qnorm(1 - pvals)
  2 * pnorm(-abs(sum(z) / sqrt(length(z))))
}

cct_combine <- function(pvals) {
  pvals <- pmax(pvals, P_CLIP)
  t_stat <- mean(tan((0.5 - pvals) * pi))
  0.5 - atan(t_stat) / pi
}

# ── Consensus Methods ────────────────────────────────────────────────────────
run_consensus <- function(results, direction = "up") {
  genes    <- results[[1]]$gene
  n_genes  <- length(genes)
  n_methods <- length(results)

  rank_lists <- vector("list", n_methods)
  for (i in seq_len(n_methods)) {
    p_one_val <- p_one(results[[i]]$p, results[[i]]$logFC, direction)
    rank_lists[[i]] <- results[[i]]$gene[order(p_one_val)]
  }

  rra_out <- tryCatch(
    aggregateRanks(rank_lists, N = n_genes, method = "RRA"),
    error = function(e) NULL
  )
  if (is.null(rra_out)) rra_out <- data.frame(Name = genes, Score = rep(1, n_genes))

  col_name <- function(suffix) paste0(direction, "_", suffix)

  tibble(
    gene = rra_out$Name,
    !!col_name("rra_p")   := rra_out$Score,
    !!col_name("rra_fdr") := p.adjust(rra_out$Score, "BH"),
    !!col_name("cct_p")   := apply(do.call(cbind, lapply(results, function(x) {
      p_one(x$p, x$logFC, direction)
    })), 1, cct_combine),
    !!col_name("cct_fdr") := p.adjust(.data[[col_name("cct_p")]], "BH"),
    !!col_name("fisher_p") := apply(do.call(cbind, lapply(results, function(x) {
      p_one(x$p, x$logFC, direction)
    })), 1, fisher_combine),
    !!col_name("fisher_fdr") := p.adjust(.data[[col_name("fisher_p")]], "BH"),
    !!col_name("stouffer_p") := apply(do.call(cbind, lapply(results, function(x) {
      p_one(x$p, x$logFC, direction)
    })), 1, stouffer_combine),
    !!col_name("stouffer_fdr") := p.adjust(.data[[col_name("stouffer_p")]], "BH")
  )
}

# ── Evaluation helpers ───────────────────────────────────────────────────────

compute_metrics <- function(pred_df, truth_vec, method_label) {
  # pred_df: columns gene, significant (logical), direction (integer)
  gene_names <- pred_df$gene
  truth_sub  <- truth_vec[gene_names]
  is_sig     <- pred_df$significant & !is.na(pred_df$significant)
  pred_dir   <- pred_df$direction
  pred_dir[is.na(pred_dir)] <- 0L

  tp <- sum(is_sig & (truth_sub != 0) & (pred_dir == truth_sub), na.rm = TRUE)
  fp <- sum(is_sig & ((truth_sub == 0) | (pred_dir != truth_sub)), na.rm = TRUE)
  fn <- sum((!is_sig) & (truth_sub != 0), na.rm = TRUE)
  n_de_true <- sum(truth_sub != 0, na.rm = TRUE)

  precision <- if (tp + fp > 0) tp / (tp + fp) else NA_real_
  recall    <- if (n_de_true > 0) tp / n_de_true else NA_real_
  f1        <- if (!is.na(precision) && !is.na(recall) && (precision + recall) > 0) {
    2 * precision * recall / (precision + recall)
  } else NA_real_

  data.frame(method = method_label, tp = tp, fp = fp, fn = fn,
             n_de_true = n_de_true, precision = precision, recall = recall,
             f1 = f1, stringsAsFactors = FALSE)
}

compute_pr_points <- function(scores_df, score_col, truth_vec, method_label,
                              direction_col = NULL) {
  gene_names <- scores_df$gene
  truth_sub  <- truth_vec[gene_names]
  n_de       <- sum(truth_sub != 0, na.rm = TRUE)

  ord       <- order(scores_df[[score_col]], na.last = TRUE)
  ord_genes <- gene_names[ord]

  ranks_to_eval <- c(10L, 20L, 50L, 100L, 200L, 500L, 1000L, 2000L, 5000L, N_GENES)
  ranks_to_eval <- ranks_to_eval[ranks_to_eval <= N_GENES]

  out <- data.frame()
  for (k in ranks_to_eval) {
    top_genes  <- ord_genes[1:k]
    top_truth  <- truth_sub[top_genes]
    if (!is.null(direction_col)) {
      top_dirs <- scores_df[[direction_col]][ord[1:k]]
      tp <- sum((top_truth != 0) & (top_dirs == top_truth), na.rm = TRUE)
    } else {
      tp <- sum(top_truth != 0, na.rm = TRUE)
    }
    fp <- k - tp
    precision <- if (k > 0) tp / k else NA_real_
    recall    <- if (n_de > 0) tp / n_de else NA_real_
    out <- rbind(out, data.frame(
      method = method_label, rank_k = k,
      precision = precision, recall = recall, stringsAsFactors = FALSE
    ))
  }
  out
}

# ═══════════════════════════════════════════════════════════════════════════════
# 1. NULL CALIBRATION  (500 repeats, 0% DEG)
# ═══════════════════════════════════════════════════════════════════════════════
cat("=== Null Simulation (0% DEG, 500 repeats) ===\n")
null_results <- data.frame()
t0 <- Sys.time()

for (sim in 1:500) {
  sim_data <- simulate_counts(deg_pct = 0)

  deseq2_res <- tryCatch(run_deseq2(sim_data$counts, N_REP), error = function(e) NULL)
  edger_res  <- tryCatch(run_edger(sim_data$counts, N_REP),  error = function(e) NULL)
  limma_res  <- tryCatch(run_limma(sim_data$counts, N_REP),  error = function(e) NULL)

  if (any(sapply(list(deseq2_res, edger_res, limma_res), is.null))) next

  res_list <- list(DESeq2 = deseq2_res, edgeR = edger_res, limma = limma_res)

  # Single methods
  for (mname in names(res_list)) {
    n_sig <- sum(res_list[[mname]]$padj < ALPHA, na.rm = TRUE)
    null_results <- rbind(null_results, data.frame(
      sim_id = sim, method = mname, n_significant = n_sig,
      type_I_error = n_sig / N_GENES, stringsAsFactors = FALSE
    ))
  }

  # Consensus
  cr_up   <- run_consensus(res_list, "up")
  cr_down <- run_consensus(res_list, "down")

  cons <- cr_up %>%
    inner_join(cr_down %>% select(gene, down_rra_fdr, down_cct_fdr,
                                  down_fisher_fdr, down_stouffer_fdr), by = "gene") %>%
    mutate(
      best_rra_fdr      = pmin(up_rra_fdr, down_rra_fdr, na.rm = TRUE),
      best_cct_fdr      = pmin(up_cct_fdr, down_cct_fdr, na.rm = TRUE),
      best_fisher_fdr   = pmin(up_fisher_fdr, down_fisher_fdr, na.rm = TRUE),
      best_stouffer_fdr = pmin(up_stouffer_fdr, down_stouffer_fdr, na.rm = TRUE)
    )

  # RRA
  n_sig <- sum(cons$best_rra_fdr < ALPHA, na.rm = TRUE)
  null_results <- rbind(null_results, data.frame(
    sim_id = sim, method = "RRA", n_significant = n_sig,
    type_I_error = n_sig / N_GENES, stringsAsFactors = FALSE))

  # CCT
  n_sig <- sum(cons$best_cct_fdr < ALPHA, na.rm = TRUE)
  null_results <- rbind(null_results, data.frame(
    sim_id = sim, method = "CCT", n_significant = n_sig,
    type_I_error = n_sig / N_GENES, stringsAsFactors = FALSE))

  # Fisher
  n_sig <- sum(cons$best_fisher_fdr < ALPHA, na.rm = TRUE)
  null_results <- rbind(null_results, data.frame(
    sim_id = sim, method = "Fisher", n_significant = n_sig,
    type_I_error = n_sig / N_GENES, stringsAsFactors = FALSE))

  # Stouffer
  n_sig <- sum(cons$best_stouffer_fdr < ALPHA, na.rm = TRUE)
  null_results <- rbind(null_results, data.frame(
    sim_id = sim, method = "Stouffer", n_significant = n_sig,
    type_I_error = n_sig / N_GENES, stringsAsFactors = FALSE))

  # CCT-RRA (Tier A-like)
  cct_rra_sig <- (cons$up_rra_fdr < ALPHA & cons$up_cct_fdr < ALPHA) |
                 (cons$down_rra_fdr < ALPHA & cons$down_cct_fdr < ALPHA)
  n_sig <- sum(cct_rra_sig, na.rm = TRUE)
  null_results <- rbind(null_results, data.frame(
    sim_id = sim, method = "CCT-RRA", n_significant = n_sig,
    type_I_error = n_sig / N_GENES, stringsAsFactors = FALSE))

  if (sim %% 100 == 0) {
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    cat(sprintf("  Sim %d/500 (%.0fs elapsed)\n", sim, elapsed))
  }
}

# Summarise null
null_summary <- null_results %>%
  group_by(method) %>%
  summarise(
    n_sims        = n(),
    mean_type_I   = mean(type_I_error, na.rm = TRUE),
    se_type_I     = sd(type_I_error, na.rm = TRUE) / sqrt(n()),
    ci95_low      = mean_type_I - 1.96 * se_type_I,
    ci95_high     = mean_type_I + 1.96 * se_type_I,
    median_type_I = median(type_I_error, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nNull calibration summary:\n")
print(as.data.frame(null_summary), digits = 4, row.names = FALSE)

write_tsv(null_results, file.path(OUT_DIR, "null_calibration.tsv"))
cat(sprintf("Saved: null_calibration.tsv  (%d rows)\n", nrow(null_results)))

# ═══════════════════════════════════════════════════════════════════════════════
# 2. NON-NULL SIMULATION  (4 scenarios × 50 repeats)
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== Non-Null Simulation ===\n")

scenarios <- list(
  list(deg_pct = 0.05, log2fc = 0.5, label = "5%_FC0.5"),
  list(deg_pct = 0.05, log2fc = 1.0, label = "5%_FC1.0"),
  list(deg_pct = 0.10, log2fc = 1.0, label = "10%_FC1.0"),
  list(deg_pct = 0.20, log2fc = 1.0, label = "20%_FC1.0")
)

set.seed(123)

nonnull_results <- data.frame()
pr_curve_data   <- data.frame()

for (sc in scenarios) {
  cat(sprintf("\n--- %s (DEG=%.0f%%, log2FC=%.1f) ---\n",
              sc$label, sc$deg_pct * 100, sc$log2fc))
  t1 <- Sys.time()

  for (rep in 1:50) {
    sim_data  <- simulate_counts(deg_pct = sc$deg_pct, log2fc = sc$log2fc)
    truth_vec <- setNames(sim_data$truth, rownames(sim_data$counts))

    deseq2_res <- tryCatch(run_deseq2(sim_data$counts, N_REP), error = function(e) NULL)
    edger_res  <- tryCatch(run_edger(sim_data$counts, N_REP),  error = function(e) NULL)
    limma_res  <- tryCatch(run_limma(sim_data$counts, N_REP),  error = function(e) NULL)

    if (any(sapply(list(deseq2_res, edger_res, limma_res), is.null))) next

    res_list <- list(DESeq2 = deseq2_res, edgeR = edger_res, limma = limma_res)

    # ── Single methods ──
    for (mname in names(res_list)) {
      df <- res_list[[mname]]
      df$significant <- df$padj < ALPHA & !is.na(df$padj)
      df$direction   <- sign(df$logFC)
      df$direction[is.na(df$direction)] <- 0L

      met <- compute_metrics(df, truth_vec, mname)
      met$sim_id  <- rep; met$deg_pct <- sc$deg_pct
      met$log2fc  <- sc$log2fc; met$scenario <- sc$label
      nonnull_results <- rbind(nonnull_results, met)

      pr <- compute_pr_points(df, "p", truth_vec, mname, "direction")
      pr$sim_id  <- rep; pr$deg_pct <- sc$deg_pct
      pr$log2fc  <- sc$log2fc; pr$scenario <- sc$label
      pr_curve_data <- rbind(pr_curve_data, pr)
    }

    # ── Consensus methods ──
    cr_up   <- run_consensus(res_list, "up")
    cr_down <- run_consensus(res_list, "down")

    cons <- cr_up %>%
      inner_join(cr_down %>% select(gene, down_rra_fdr, down_cct_fdr,
                                    down_fisher_fdr, down_stouffer_fdr), by = "gene") %>%
      mutate(
        best_rra_fdr      = pmin(up_rra_fdr, down_rra_fdr, na.rm = TRUE),
        best_cct_fdr      = pmin(up_cct_fdr, down_cct_fdr, na.rm = TRUE),
        best_fisher_fdr   = pmin(up_fisher_fdr, down_fisher_fdr, na.rm = TRUE),
        best_stouffer_fdr = pmin(up_stouffer_fdr, down_stouffer_fdr, na.rm = TRUE),
        rra_dir           = ifelse(up_rra_fdr < down_rra_fdr, 1L, -1L),
        cct_dir           = ifelse(up_cct_fdr < down_cct_fdr, 1L, -1L),
        fisher_dir        = ifelse(up_fisher_fdr < down_fisher_fdr, 1L, -1L),
        stouffer_dir      = ifelse(up_stouffer_fdr < down_stouffer_fdr, 1L, -1L)
      )

    # ── RRA ──
    df <- data.frame(gene = cons$gene,
                     significant = cons$best_rra_fdr < ALPHA,
                     direction = cons$rra_dir, stringsAsFactors = FALSE)
    met_row <- compute_metrics(df, truth_vec, "RRA")
    met_row$sim_id <- rep; met_row$deg_pct <- sc$deg_pct
    met_row$log2fc <- sc$log2fc; met_row$scenario <- sc$label
    tmp <- cons %>% mutate(rank_score = best_rra_fdr, pr_dir = rra_dir)
    pr_row <- compute_pr_points(tmp, "rank_score", truth_vec, "RRA", "pr_dir")
    pr_row$sim_id <- rep; pr_row$deg_pct <- sc$deg_pct
    pr_row$log2fc <- sc$log2fc; pr_row$scenario <- sc$label
    nonnull_results <- rbind(nonnull_results, met_row)
    pr_curve_data   <- rbind(pr_curve_data, pr_row)

    # ── CCT ──
    df <- data.frame(gene = cons$gene,
                     significant = cons$best_cct_fdr < ALPHA,
                     direction = cons$cct_dir, stringsAsFactors = FALSE)
    met_row <- compute_metrics(df, truth_vec, "CCT")
    met_row$sim_id <- rep; met_row$deg_pct <- sc$deg_pct
    met_row$log2fc <- sc$log2fc; met_row$scenario <- sc$label
    tmp <- cons %>% mutate(rank_score = best_cct_fdr, pr_dir = cct_dir)
    pr_row <- compute_pr_points(tmp, "rank_score", truth_vec, "CCT", "pr_dir")
    pr_row$sim_id <- rep; pr_row$deg_pct <- sc$deg_pct
    pr_row$log2fc <- sc$log2fc; pr_row$scenario <- sc$label
    nonnull_results <- rbind(nonnull_results, met_row)
    pr_curve_data   <- rbind(pr_curve_data, pr_row)

    # ── Fisher ──
    df <- data.frame(gene = cons$gene,
                     significant = cons$best_fisher_fdr < ALPHA,
                     direction = cons$fisher_dir, stringsAsFactors = FALSE)
    met_row <- compute_metrics(df, truth_vec, "Fisher")
    met_row$sim_id <- rep; met_row$deg_pct <- sc$deg_pct
    met_row$log2fc <- sc$log2fc; met_row$scenario <- sc$label
    tmp <- cons %>% mutate(rank_score = best_fisher_fdr, pr_dir = fisher_dir)
    pr_row <- compute_pr_points(tmp, "rank_score", truth_vec, "Fisher", "pr_dir")
    pr_row$sim_id <- rep; pr_row$deg_pct <- sc$deg_pct
    pr_row$log2fc <- sc$log2fc; pr_row$scenario <- sc$label
    nonnull_results <- rbind(nonnull_results, met_row)
    pr_curve_data   <- rbind(pr_curve_data, pr_row)

    # ── Stouffer ──
    df <- data.frame(gene = cons$gene,
                     significant = cons$best_stouffer_fdr < ALPHA,
                     direction = cons$stouffer_dir, stringsAsFactors = FALSE)
    met_row <- compute_metrics(df, truth_vec, "Stouffer")
    met_row$sim_id <- rep; met_row$deg_pct <- sc$deg_pct
    met_row$log2fc <- sc$log2fc; met_row$scenario <- sc$label
    tmp <- cons %>% mutate(rank_score = best_stouffer_fdr, pr_dir = stouffer_dir)
    pr_row <- compute_pr_points(tmp, "rank_score", truth_vec, "Stouffer", "pr_dir")
    pr_row$sim_id <- rep; pr_row$deg_pct <- sc$deg_pct
    pr_row$log2fc <- sc$log2fc; pr_row$scenario <- sc$label
    nonnull_results <- rbind(nonnull_results, met_row)
    pr_curve_data   <- rbind(pr_curve_data, pr_row)

    # CCT-RRA (Tier A)
    cons <- cons %>% mutate(
      cct_rra_up   = up_rra_fdr < ALPHA & up_cct_fdr < ALPHA,
      cct_rra_down = down_rra_fdr < ALPHA & down_cct_fdr < ALPHA,
      tier_sig     = cct_rra_up | cct_rra_down,
      tier_dir     = case_when(
        cct_rra_up & !cct_rra_down ~ 1L,
        !cct_rra_up & cct_rra_down ~ -1L,
        cct_rra_up & cct_rra_down ~ ifelse(up_rra_fdr < down_rra_fdr, 1L, -1L),
        TRUE ~ 0L
      ),
      tier_score = pmin(
        ifelse(cct_rra_up,   pmax(up_rra_fdr, up_cct_fdr, na.rm = TRUE), NA_real_),
        ifelse(cct_rra_down, pmax(down_rra_fdr, down_cct_fdr, na.rm = TRUE), NA_real_),
        na.rm = TRUE
      )
    )

    df_tier <- data.frame(gene = cons$gene, significant = cons$tier_sig,
                          direction = cons$tier_dir, stringsAsFactors = FALSE)
    met <- compute_metrics(df_tier, truth_vec, "CCT-RRA")
    met$sim_id <- rep; met$deg_pct <- sc$deg_pct
    met$log2fc <- sc$log2fc; met$scenario <- sc$label
    nonnull_results <- rbind(nonnull_results, met)

    tmp <- cons %>% mutate(rank_score = tier_score)
    pr <- compute_pr_points(tmp, "rank_score", truth_vec, "CCT-RRA", "tier_dir")
    pr$sim_id <- rep; pr$deg_pct <- sc$deg_pct
    pr$log2fc <- sc$log2fc; pr$scenario <- sc$label
    pr_curve_data <- rbind(pr_curve_data, pr)

    if (rep %% 10 == 0) cat(sprintf("  %d/50\n", rep))
  }
  elapsed <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  cat(sprintf("  Done in %.0fs\n", elapsed))
}

# Save nonnull
write_tsv(nonnull_results, file.path(OUT_DIR, "nonnull_results.tsv"))
write_tsv(pr_curve_data,   file.path(OUT_DIR, "pr_curve_data.tsv"))
cat(sprintf("\nSaved: nonnull_results.tsv (%d rows)\n", nrow(nonnull_results)))
cat(sprintf("Saved: pr_curve_data.tsv   (%d rows)\n", nrow(pr_curve_data)))

# ═══════════════════════════════════════════════════════════════════════════════
# 3. TABLE 2 – Aggregated manuscript data
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== Table 2 – Aggregated Results ===\n")

table2_data <- nonnull_results %>%
  group_by(scenario, deg_pct, log2fc, method) %>%
  summarise(
    n_reps         = n(),
    mean_precision = mean(precision, na.rm = TRUE),
    sd_precision   = sd(precision, na.rm = TRUE),
    mean_recall    = mean(recall, na.rm = TRUE),
    sd_recall      = sd(recall, na.rm = TRUE),
    mean_f1        = mean(f1, na.rm = TRUE),
    sd_f1          = sd(f1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(scenario, desc(mean_f1))

cat("\nPerformance by scenario and method:\n")
print(as.data.frame(table2_data), digits = 4, row.names = FALSE)

write_tsv(table2_data, file.path(OUT_DIR, "table2_data.tsv"))
cat(sprintf("Saved: table2_data.tsv (%d rows)\n", nrow(table2_data)))

# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== Benchmark Complete ===\n")
cat(sprintf("Null sims:          %d (target 500)\n",
    max(null_results$sim_id, 0)))
cat(sprintf("Non-null sims:      %d (target 200)\n",
    nrow(distinct(nonnull_results, sim_id, scenario))))
cat(sprintf("Output directory:   %s\n", OUT_DIR))
cat(sprintf("Null calibration:   %s\n", file.path(OUT_DIR, "null_calibration.tsv")))
cat(sprintf("Nonnull results:    %s\n", file.path(OUT_DIR, "nonnull_results.tsv")))
cat(sprintf("PR curve data:      %s\n", file.path(OUT_DIR, "pr_curve_data.tsv")))
cat(sprintf("Table 2 data:       %s\n", file.path(OUT_DIR, "table2_data.tsv")))
