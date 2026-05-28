#!/usr/bin/env Rscript
# workflow/scripts/benchmark_subsampling.R
# Exhaustive 9-subset stability benchmark for multi-quantifier consensus DEA
# Design: Drosophila 60d_vs_1d, 3 reps/group, choose 2 per group = 9 subsets
# Standalone script — replicates consensus functions from run_consensus_dea.R

# ── Packages ──────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(DESeq2)
  library(RobustRankAggreg)
  library(ggplot2)
  library(matrixStats)
})

# ── Config ────────────────────────────────────────────────────────────────────
FDR_THRESHOLD      <- 0.05
LFC_THRESHOLD      <- 1.0
P_CLIP             <- 1e-16
CONTRAST_TREAT     <- "Wolbachia_infected"
CONTRAST_CTRL      <- "Wolbachia_free"
QUANTIFIERS        <- c("featurecounts", "salmon", "stringtie", "kallisto")
TOP_N_VALUES       <- c(100, 250, 500, 1000)
FIXED_DIRECTION    <- "up"  # stability calculated within up-regulated genes

# Tier thresholds from config.yaml — CCT max_fdr set to 1.0 so consensus is
# RRA-driven (per plan: "Set CCT threshold to 1.0 for CCT")
TIERS <- list(
  tier_a = list(min_support = 4L, min_sign_consistency = 4L,
                max_rra_fdr = 0.05, max_cct_fdr = 1.0, max_logfc_cv = 1.0),
  tier_b = list(min_support = 3L, min_sign_consistency = 3L,
                max_rra_fdr = 0.10, max_cct_fdr = 1.0, max_logfc_cv = 1.25),
  tier_c = list(min_support = 2L, min_sign_consistency = 2L,
                max_rra_fdr = 0.25, max_cct_fdr = 1.0, max_logfc_cv = 1.50)
)

# ── File paths ────────────────────────────────────────────────────────────────
MATRIX_DIR        <- "results/05.quantification/matrices"
FC_MATRIX_PATH    <- file.path(MATRIX_DIR, "featurecounts",
                                "featurecounts_gene_counts_matrix.tsv")
SALMON_PATH       <- file.path(MATRIX_DIR, "salmon",
                                "salmon_gene_counts_matrix.tsv")
STRINGTIE_PATH    <- file.path(MATRIX_DIR, "stringtie",
                                "stringtie_gene_counts_matrix.tsv")
KALLISTO_PATH     <- file.path(MATRIX_DIR, "kallisto",
                                "kallisto_transcript_counts_matrix.tsv")
TX2GENE_PATH      <- "results/00.reference/tx2gene_master.tsv"
SAMPLES_PATH      <- "data/fastq/samples.tsv"
OUT_DIR           <- "results/benchmark"
FIG_DIR           <- file.path(OUT_DIR, "figures")

# ══════════════════════════════════════════════════════════════════════════════
# REPLICATED CONSENSUS FUNCTIONS — from run_consensus_dea.R
# ══════════════════════════════════════════════════════════════════════════════

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

as_flag <- function(x, default = FALSE) {
  if (is.null(x) || length(x) == 0 || is.na(x[1])) return(default)
  if (is.logical(x)) return(x)
  tolower(as.character(x)) %in% c("true", "1", "yes")
}

clip_probabilities <- function(p, eps) {
  pmin(pmax(p, eps), 1 - eps)
}

build_rank_list <- function(df, direction) {
  ranked <- df %>%
    filter(!is.na(.data$gene_id_standard), !is.na(.data$P.Value),
           !is.na(.data$logFC))
  if (direction == "up") {
    ranked <- ranked %>%
      filter(.data$logFC > 0) %>%
      arrange(.data$P.Value, desc(.data$logFC), .data$gene_id_standard)
  } else {
    ranked <- ranked %>%
      filter(.data$logFC < 0) %>%
      arrange(.data$P.Value, .data$logFC, .data$gene_id_standard)
  }
  unique(ranked$gene_id_standard)
}

compute_rra_scores <- function(data_list, universe, direction) {
  rank_lists <- lapply(data_list, build_rank_list, direction = direction)
  rank_lists <- rank_lists[lengths(rank_lists) > 0]

  result <- tibble(gene_id_standard = universe, p = rep(1, length(universe)))
  if (length(rank_lists) == 0) {
    result$fdr <- p.adjust(result$p, method = "BH")
    return(result)
  }

  rra_out <- RobustRankAggreg::aggregateRanks(
    glist = rank_lists, N = length(universe), method = "RRA", exact = FALSE
  )
  observed <- tibble(gene_id_standard = rra_out$Name, p_observed = rra_out$Score)
  result <- result %>%
    left_join(observed, by = "gene_id_standard") %>%
    mutate(
      p = coalesce(.data$p_observed, .data$p),
      p = clip_probabilities(.data$p, eps = 1e-300)
    ) %>%
    select(-any_of("p_observed"))
  result$fdr <- p.adjust(result$p, method = "BH")
  result
}

compute_cct_scores <- function(logfc_matrix, p_matrix, universe, direction, eps) {
  mask <- !is.na(logfc_matrix) & !is.na(p_matrix)
  if (direction == "up") {
    mask <- mask & (logfc_matrix > 0)
  } else {
    mask <- mask & (logfc_matrix < 0)
  }

  valid_p <- p_matrix
  valid_p[!mask] <- NA
  valid_p <- clip_probabilities(valid_p, eps = eps)

  tan_mat <- tan((0.5 - valid_p) * pi)
  t_stat <- rowMeans(tan_mat, na.rm = TRUE)
  cct_p <- 0.5 - atan(t_stat) / pi
  cct_p[is.nan(cct_p)] <- 1
  cct_p <- pmin(pmax(cct_p, eps), 1)

  tibble(gene_id_standard = universe, p = cct_p,
         fdr = p.adjust(cct_p, method = "BH"))
}

choose_direction <- function(up_support, down_support) {
  if (up_support == 0 && down_support == 0) return("none")
  if (up_support == down_support) return("mixed")
  if (up_support > 0 && down_support > 0) {
    if (up_support >= 2 * down_support) return("up")
    if (down_support >= 2 * up_support) return("down")
    return("mixed")
  }
  if (up_support > down_support) return("up")
  if (down_support > up_support) return("down")
  "none"
}

compute_consensus_logfc <- function(logfc_matrix, consensus_direction) {
  res <- numeric(nrow(logfc_matrix))
  for (i in seq_len(nrow(logfc_matrix))) {
    values <- logfc_matrix[i, ]
    dir_i <- consensus_direction[i]
    if (dir_i == "mixed" || dir_i == "none") { res[i] <- NA_real_; next }
    values <- values[!is.na(values)]
    if (length(values) == 0) { res[i] <- NA_real_; next }
    if (dir_i == "up") values <- values[values > 0]
    else if (dir_i == "down") values <- values[values < 0]
    if (length(values) == 0) res[i] <- NA_real_
    else res[i] <- median(values)
  }
  res
}

compute_logfc_cv <- function(logfc_matrix, consensus_direction) {
  res <- numeric(nrow(logfc_matrix))
  for (i in seq_len(nrow(logfc_matrix))) {
    values <- logfc_matrix[i, ]
    dir_i <- consensus_direction[i]
    if (dir_i == "mixed" || dir_i == "none") { res[i] <- NA_real_; next }
    values <- values[!is.na(values)]
    if (length(values) < 2) { res[i] <- NA_real_; next }
    if (dir_i == "up") values <- values[values > 0]
    else if (dir_i == "down") values <- abs(values[values < 0])
    else values <- abs(values)
    if (length(values) < 2) { res[i] <- NA_real_; next }
    mu <- mean(values)
    if (is.na(mu) || mu == 0) { res[i] <- NA_real_ }
    else res[i] <- sd(values) / mu
  }
  res
}

assign_tier <- function(support_n, sign_consistency_n, best_rra_fdr,
                        best_cct_fdr, logfc_cv, tiers, consensus_direction = "up") {
  if (consensus_direction == "mixed" || consensus_direction == "none")
    return("unclassified")
  if (support_n < 2) return("unclassified")

  ta <- tiers$tier_a; tb <- tiers$tier_b; tc <- tiers$tier_c

  if (support_n >= ta$min_support &&
      sign_consistency_n >= ta$min_sign_consistency &&
      sign_consistency_n == support_n &&
      best_rra_fdr <= ta$max_rra_fdr &&
      best_cct_fdr <= ta$max_cct_fdr &&
      (is.na(logfc_cv) || logfc_cv <= ta$max_logfc_cv))
    return("Tier_A")

  if (support_n >= tb$min_support &&
      sign_consistency_n >= tb$min_sign_consistency &&
      best_rra_fdr <= tb$max_rra_fdr &&
      best_cct_fdr <= tb$max_cct_fdr &&
      (is.na(logfc_cv) || logfc_cv <= tb$max_logfc_cv))
    return("Tier_B")

  if (support_n >= tc$min_support &&
      sign_consistency_n >= tc$min_sign_consistency &&
      best_rra_fdr <= tc$max_rra_fdr &&
      best_cct_fdr <= tc$max_cct_fdr &&
      (is.na(logfc_cv) || logfc_cv <= tc$max_logfc_cv))
    return("Tier_C")

  "unclassified"
}

# ══════════════════════════════════════════════════════════════════════════════
# DATA LOADING
# ══════════════════════════════════════════════════════════════════════════════

#' Load a gene-level count matrix (featureCounts, Salmon, or StringTie format)
#' Returns: tibble with gene_id column + sample columns, values coerced to integer
load_gene_matrix <- function(path, gene_col) {
  cat(sprintf("  Loading: %s\n", basename(path)))
  mat <- read_tsv(path, show_col_types = FALSE)
  # Rename first column to standard name
  names(mat)[1] <- "gene_id"
  # Drop any other non-sample columns (e.g., chr, start, end, strand, length)
  sample_cols <- intersect(colnames(mat)[-1],
                           c("sub_SRR14101759","sub_SRR14101760","sub_SRR14101761",
                             "sub_SRR14101762","sub_SRR14101763","sub_SRR14101764"))
  if (length(sample_cols) == 0) {
    # Fallback: all non-gene_id columns
    sample_cols <- setdiff(colnames(mat), "gene_id")
  }
  mat <- mat %>% select("gene_id", all_of(sample_cols))
  # Convert count columns to numeric
  mat <- mat %>%
    mutate(across(all_of(sample_cols), ~ coalesce(as.numeric(.x), 0)))
  mat
}

#' Load kallisto transcript-level matrix, map to gene-level via tx2gene_master
load_kallisto_as_gene <- function(kallisto_path, tx2gene_path) {
  cat(sprintf("  Loading: %s (transcript → gene mapping)\n",
              basename(kallisto_path)))
  kt <- read_tsv(kallisto_path, show_col_types = FALSE)
  # Drop length/eff_length metadata columns
  meta_cols <- intersect(colnames(kt), c("length", "eff_length"))
  sample_cols <- setdiff(colnames(kt), c("target_id", meta_cols))
  kt <- kt %>% select("target_id", all_of(sample_cols))
  kt <- kt %>%
    mutate(across(all_of(sample_cols), ~ coalesce(as.numeric(.x), 0)))

  # Load tx2gene mapping
  tx2g <- read_tsv(tx2gene_path, show_col_types = FALSE)
  # Filter to kallisto-appropriate entries (use "reference" entries with resolved IDs)
  tx2g <- tx2g %>%
    filter(.data$quantifier %in% c("reference"),
           as_flag(.data$allow_consensus_main)) %>%
    mutate(gene_id = coalesce(.data$gene_id_resolved, .data$gene_id_original)) %>%
    filter(!is.na(.data$gene_id), .data$gene_id != "") %>%
    select("transcript_id", "gene_id")

  # Join and sum by gene
  kt_mapped <- kt %>%
    left_join(tx2g, by = c("target_id" = "transcript_id")) %>%
    filter(!is.na(.data$gene_id)) %>%
    select(-"target_id")

  # Sum counts per gene
  kt_mapped <- kt_mapped %>%
    group_by(.data$gene_id) %>%
    summarise(across(everything(), sum), .groups = "drop")

  kt_mapped
}

# ══════════════════════════════════════════════════════════════════════════════
# SUBSET UTILITIES
# ══════════════════════════════════════════════════════════════════════════════

#' Generate exhaustive subsets: choose 2 from each group
#' Returns list of data frames, each with 4 samples (2 treated, 2 control)
generate_subsets <- function(samples_df) {
  treated <- samples_df %>% filter(.data$group == CONTRAST_TREAT) %>% pull("sample")
  control <- samples_df %>% filter(.data$group == CONTRAST_CTRL) %>% pull("sample")

  treated_pairs <- combn(treated, 2, simplify = FALSE)
  control_pairs <- combn(control, 2, simplify = FALSE)

  subsets <- list()
  idx <- 1
  for (tp in treated_pairs) {
    for (cp in control_pairs) {
      subsets[[idx]] <- c(tp, cp)
      idx <- idx + 1
    }
  }
  subsets
}

#' Run DESeq2 on a subset of samples for one quantifier's count matrix
#' Returns: tibble with gene_id, gene_id_standard, logFC, P.Value, adj.P.Val
run_deseq2_subset <- function(count_mat, subset_samples, quantifier_name) {
  # Subset to selected samples
  available_samps <- intersect(subset_samples, colnames(count_mat))
  if (length(available_samps) < 4) {
    warning(sprintf("  [%s] Only %d samples available, skipping",
                    quantifier_name, length(available_samps)))
    return(tibble(gene_id_standard = character(), gene_id = character(),
                  logFC = double(), P.Value = double(), adj.P.Val = double()))
  }

  cnt <- count_mat %>%
    select("gene_id", all_of(available_samps)) %>%
    tibble::column_to_rownames("gene_id")

  # Remove genes with all-zero counts in this subset
  cnt <- cnt[rowSums(cnt) > 0, , drop = FALSE]
  if (nrow(cnt) == 0) {
    return(tibble(gene_id_standard = character(), gene_id = character(),
                  logFC = double(), P.Value = double(), adj.P.Val = double()))
  }

  # Round to integers (Salmon gives fractional estimates)
  cnt <- round(as.matrix(cnt))
  storage.mode(cnt) <- "integer"

  # Design matrix
  group_labels <- ifelse(colnames(cnt) %in%
    samples_df$sample[samples_df$group == CONTRAST_TREAT],
    CONTRAST_TREAT, CONTRAST_CTRL)
  coldata <- data.frame(
    row.names = colnames(cnt),
    group = factor(group_labels, levels = c(CONTRAST_CTRL, CONTRAST_TREAT))
  )

  # DESeq2
  suppressWarnings({
    dds <- tryCatch(
      DESeqDataSetFromMatrix(countData = cnt, colData = coldata, design = ~group),
      error = function(e) NULL
    )
  })
  if (is.null(dds)) {
    return(tibble(gene_id_standard = character(), gene_id = character(),
                  logFC = double(), P.Value = double(), adj.P.Val = double()))
  }

  dds <- tryCatch(
    DESeq(dds, quiet = TRUE),
    error = function(e) NULL
  )
  if (is.null(dds)) {
    return(tibble(gene_id_standard = character(), gene_id = character(),
                  logFC = double(), P.Value = double(), adj.P.Val = double()))
  }

  res <- tryCatch(
    results(dds, contrast = c("group", CONTRAST_TREAT, CONTRAST_CTRL),
            alpha = FDR_THRESHOLD, independentFiltering = TRUE),
    error = function(e) NULL
  )
  if (is.null(res)) {
    return(tibble(gene_id_standard = character(), gene_id = character(),
                  logFC = double(), P.Value = double(), adj.P.Val = double()))
  }

  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene_id") %>%
    as_tibble() %>%
    mutate(
      quantifier = quantifier_name,
      gene_id_standard = .data$gene_id,
      logFC = .data$log2FoldChange,
      P.Value = .data$pvalue,
      adj.P.Val = .data$padj
    ) %>%
    select("gene_id_standard", "gene_id", "logFC", "P.Value", "adj.P.Val",
           "quantifier") %>%
    filter(!is.na(.data$gene_id_standard))

  res_df
}

# ══════════════════════════════════════════════════════════════════════════════
# CONSENSUS PER SUBSET
# ══════════════════════════════════════════════════════════════════════════════

#' Run full consensus pipeline on one subset's DEA results
#' Returns: list(consensus_df, voting_calls, single_calls, directional_splits)
run_consensus_subset <- function(data_list, universe) {
  n_quant <- length(data_list)
  qnames <- names(data_list)

  # Build matrices
  consensus_df <- tibble(gene_id_standard = universe)

  logfc_list <- list()
  p_list <- list()
  adj_list <- list()

  for (q in qnames) {
    df <- data_list[[q]] %>%
      transmute(
        gene_id_standard = .data$gene_id_standard,
        !!paste0("logFC__", q) := .data$logFC,
        !!paste0("P.Value__", q) := .data$P.Value,
        !!paste0("adj.P.Val__", q) := .data$adj.P.Val
      )
    consensus_df <- consensus_df %>% left_join(df, by = "gene_id_standard")
    logfc_list[[q]] <- df[[paste0("logFC__", q)]]
    p_list[[q]] <- df[[paste0("P.Value__", q)]]
    adj_list[[q]] <- df[[paste0("adj.P.Val__", q)]]
  }

  # Build numeric matrices (aligned to universe)
  logfc_matrix <- do.call(cbind, lapply(qnames, function(q) {
    consensus_df[[paste0("logFC__", q)]]
  }))
  colnames(logfc_matrix) <- qnames

  adj_matrix <- do.call(cbind, lapply(qnames, function(q) {
    consensus_df[[paste0("adj.P.Val__", q)]]
  }))
  colnames(adj_matrix) <- qnames

  p_matrix <- do.call(cbind, lapply(qnames, function(q) {
    consensus_df[[paste0("P.Value__", q)]]
  }))
  colnames(p_matrix) <- qnames

  # Significance matrix
  sig_matrix <- !is.na(logfc_matrix) & !is.na(adj_matrix) &
    (adj_matrix <= FDR_THRESHOLD) & (abs(logfc_matrix) >= LFC_THRESHOLD)

  support_n       <- rowSums(sig_matrix, na.rm = TRUE)
  up_support_n    <- rowSums(sig_matrix & logfc_matrix > 0, na.rm = TRUE)
  down_support_n  <- rowSums(sig_matrix & logfc_matrix < 0, na.rm = TRUE)
  sign_consistency_n <- pmax(up_support_n, down_support_n)

  # RRA scores
  data_for_rra <- lapply(qnames, function(q) data_list[[q]])
  names(data_for_rra) <- qnames

  rra_up   <- compute_rra_scores(data_for_rra, universe, "up")
  rra_down <- compute_rra_scores(data_for_rra, universe, "down")

  # CCT scores
  cct_up   <- compute_cct_scores(logfc_matrix, p_matrix, universe, "up", P_CLIP)
  cct_down <- compute_cct_scores(logfc_matrix, p_matrix, universe, "down", P_CLIP)

  # Merge
  # Rename before join — use setNames to avoid dplyr rename issues
  rra_up_renamed   <- setNames(rra_up,   c("gene_id_standard","rra_up_p","rra_up_fdr"))
  rra_down_renamed <- setNames(rra_down, c("gene_id_standard","rra_down_p","rra_down_fdr"))
  cct_up_renamed   <- setNames(cct_up,   c("gene_id_standard","cct_up_p","cct_up_fdr"))
  cct_down_renamed <- setNames(cct_down, c("gene_id_standard","cct_down_p","cct_down_fdr"))

  consensus_df <- consensus_df %>%
    left_join(rra_up_renamed, by = "gene_id_standard") %>%
    left_join(rra_down_renamed, by = "gene_id_standard") %>%
    left_join(cct_up_renamed, by = "gene_id_standard") %>%
    left_join(cct_down_renamed, by = "gene_id_standard")

  # Direction
  cons_dir <- vapply(seq_len(nrow(consensus_df)),
                     function(i) choose_direction(up_support_n[i], down_support_n[i]),
                     character(1))

  # best RRA/CCT FDR by direction
  best_rra_fdr <- ifelse(cons_dir == "up", consensus_df$rra_up_fdr,
                   ifelse(cons_dir == "down", consensus_df$rra_down_fdr, NA_real_))
  best_cct_fdr <- ifelse(cons_dir == "up", consensus_df$cct_up_fdr,
                   ifelse(cons_dir == "down", consensus_df$cct_down_fdr, NA_real_))

  # logFC CV
  logfc_cv <- compute_logfc_cv(logfc_matrix, cons_dir)

  # Assign tiers
  tier <- vapply(seq_len(nrow(consensus_df)), function(i) {
    assign_tier(support_n[i], sign_consistency_n[i], best_rra_fdr[i],
                best_cct_fdr[i], logfc_cv[i], TIERS, cons_dir[i])
  }, character(1))

  # ── Voting calls (>=2 quantifiers agree on direction) ──
  voting_up   <- up_support_n >= 2
  voting_down <- down_support_n >= 2
  voting_de   <- voting_up | voting_down

  # ── Single-quantifier calls ──
  single_calls <- as.data.frame(sig_matrix)
  single_calls$gene_id_standard <- universe
  single_calls <- as_tibble(single_calls)

  # ── Directional splits (2-2, 3-1, etc.) ──
  directional_splits <- tibble(
    gene_id_standard = universe,
    up_n    = up_support_n,
    down_n  = down_support_n,
    pattern = paste0(up_support_n, "up_", down_support_n, "down")
  )

  # ── Return ──
  list(
    # Consensus Tier A genes
    tier_a_genes     = consensus_df$gene_id_standard[tier == "Tier_A"],
    tier_b_genes     = consensus_df$gene_id_standard[tier == "Tier_B"],
    tier_c_genes     = consensus_df$gene_id_standard[tier == "Tier_C"],
    # Voting DE genes
    voting_de_genes  = universe[voting_de],
    voting_up_genes  = universe[voting_up],
    voting_down_genes = universe[voting_down],
    # Per-quantifier DE genes
    single_calls     = single_calls,
    # Consensus ranking (by best_rra_fdr ascending)
    consensus_rank   = consensus_df %>%
      mutate(cons_dir = cons_dir) %>%
      filter(.data$cons_dir == "up") %>%
      arrange(.data$rra_up_fdr) %>%
      pull("gene_id_standard"),
    # Voting ranking (by support_n desc, then abs consensus logFC desc)
    voting_rank = {
      median_abs_lfc <- rowMedians(abs(logfc_matrix), na.rm = TRUE)
      voting_rank_df <- consensus_df %>%
        mutate(support_n = support_n, median_abs_lfc = median_abs_lfc,
               voting_up = voting_up) %>%
        filter(voting_up) %>%
        arrange(desc(.data$support_n), desc(.data$median_abs_lfc))
      voting_rank_df$gene_id_standard
    },
    # Single-quantifier rankings (by padj ascending, up-regulated)
    single_ranks = lapply(qnames, function(q) {
      df <- data_list[[q]] %>%
        filter(!is.na(.data$adj.P.Val), !is.na(.data$logFC),
               .data$logFC > 0) %>%
        arrange(.data$adj.P.Val)
      df$gene_id_standard
    }) %>% setNames(qnames),
    # Directional splits summary
    directional_splits = directional_splits,
    # Full consensus df for logging
    n_tier_a = sum(tier == "Tier_A"),
    n_tier_b = sum(tier == "Tier_B"),
    n_tier_c = sum(tier == "Tier_C"),
    n_voting_de = sum(voting_de),
    n_voting_up = sum(voting_up),
    n_voting_down = sum(voting_down),
    n_2up_2down = sum(up_support_n == 2 & down_support_n == 2)
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# STABILITY METRICS
# ══════════════════════════════════════════════════════════════════════════════

#' Compute gene-level stability: fraction of subsets where gene is called DE
#' @param calls_matrix logical matrix: genes × subsets
#' @return numeric vector (one per gene)
compute_gene_stability <- function(calls_matrix) {
  rowMeans(calls_matrix, na.rm = TRUE)
}

#' Compute mean Jaccard index for top-N lists across all subset pairs
#' @param rank_lists list of character vectors (top-N genes per subset)
#' @return mean Jaccard across all pairs
compute_topn_jaccard <- function(rank_lists) {
  n <- length(rank_lists)
  if (n < 2) return(NA_real_)
  pairs <- combn(n, 2)
  jaccards <- vapply(seq_len(ncol(pairs)), function(p) {
    i <- pairs[1, p]; j <- pairs[2, p]
    set_i <- rank_lists[[i]]; set_j <- rank_lists[[j]]
    if (length(set_i) == 0 && length(set_j) == 0) return(1.0)
    if (length(set_i) == 0 || length(set_j) == 0) return(0.0)
    length(intersect(set_i, set_j)) / length(union(set_i, set_j))
  }, numeric(1))
  mean(jaccards)
}

#' Jaccard for a specific N with dynamic list sizing
compute_jaccard_at_n <- function(rank_lists, N) {
  top_lists <- lapply(rank_lists, function(x) head(x, N))
  compute_topn_jaccard(top_lists)
}

#' Mean number of DE genes called per subset
compute_mean_yield <- function(calls_matrix) {
  mean(colSums(calls_matrix, na.rm = TRUE))
}

# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

cat("═══════════════════════════════════════════════════════════════\n")
cat("  OmniQuant-RNA: Subsampling Stability Benchmark\n")
cat("  Design: 9 exhaustive subsets (3 choose 2 per group)\n")
cat("  Contrast:", CONTRAST_TREAT, "vs", CONTRAST_CTRL, "\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# ── Step 1: Load samples ─────────────────────────────────────────────────────
cat("[1/7] Loading sample metadata...\n")
  samples_df <- read_tsv(SAMPLES_PATH, show_col_types = FALSE,
                          col_types = cols(sample = col_character(),
                                          fq1 = col_character(),
                                          fq2 = col_character(),
                                          group = col_character()))
cat(sprintf("  %d samples across %d groups\n",
            nrow(samples_df), length(unique(samples_df$group))))

# Validate groups
stopifnot(CONTRAST_TREAT %in% samples_df$group)
stopifnot(CONTRAST_CTRL %in% samples_df$group)
cat(sprintf("  %s: %d samples | %s: %d samples\n",
  CONTRAST_TREAT, sum(samples_df$group == CONTRAST_TREAT),
  CONTRAST_CTRL, sum(samples_df$group == CONTRAST_CTRL)))

# ── Step 2: Load count matrices ──────────────────────────────────────────────
cat("\n[2/7] Loading count matrices...\n")

count_mats <- list()
count_mats[["featurecounts"]] <- load_gene_matrix(FC_MATRIX_PATH, "Geneid")
count_mats[["salmon"]]        <- load_gene_matrix(SALMON_PATH, "gene_id")
count_mats[["stringtie"]]     <- load_gene_matrix(STRINGTIE_PATH, "Gene_ID")
count_mats[["kallisto"]]      <- load_kallisto_as_gene(KALLISTO_PATH, TX2GENE_PATH)

# Report dimensions
for (q in names(count_mats)) {
  cat(sprintf("  %s: %d genes × %d samples\n",
              q, nrow(count_mats[[q]]), ncol(count_mats[[q]]) - 1))
}

# ── Step 3: Compute common gene universe ─────────────────────────────────────
cat("\n[3/7] Computing common gene universe...\n")

gene_sets <- lapply(count_mats, function(m) unique(m$gene_id))
universe <- Reduce(intersect, gene_sets)
cat(sprintf("  Per-quantifier genes: %s\n",
            paste(sapply(gene_sets, length), collapse = ", ")))
cat(sprintf("  Common universe: %d genes\n", length(universe)))

if (length(universe) == 0) {
  stop("No common genes found across quantifiers. Check count matrices.")
}

# ── Step 4: Generate subsets ─────────────────────────────────────────────────
cat("\n[4/7] Generating 9 exhaustive subsets...\n")
subsets <- generate_subsets(samples_df)
cat(sprintf("  Generated %d subsets (3 choose 2 per group)\n", length(subsets)))

# ── Step 5: Run DESeq2 + consensus for each subset ───────────────────────────
cat("\n[5/7] Running DESeq2 and consensus for each subset...\n")

# Storage for results
n_subsets <- length(subsets)

# Per-gene DE calls: one column per subset, one row per gene
genes_by_method <- list(
  tier_a         = matrix(FALSE, nrow = length(universe), ncol = n_subsets),
  voting_de      = matrix(FALSE, nrow = length(universe), ncol = n_subsets),
  featurecounts  = matrix(FALSE, nrow = length(universe), ncol = n_subsets),
  salmon         = matrix(FALSE, nrow = length(universe), ncol = n_subsets),
  stringtie      = matrix(FALSE, nrow = length(universe), ncol = n_subsets),
  kallisto       = matrix(FALSE, nrow = length(universe), ncol = n_subsets)
)
rownames(genes_by_method[["tier_a"]]) <- universe
for (m in names(genes_by_method)) {
  rownames(genes_by_method[[m]]) <- universe
}

# Rank lists for top-N Jaccard
rank_lists <- list(
  tier_a    = vector("list", n_subsets),
  voting    = vector("list", n_subsets),
  featurecounts = vector("list", n_subsets),
  salmon    = vector("list", n_subsets),
  stringtie = vector("list", n_subsets),
  kallisto  = vector("list", n_subsets)
)

# Per-subset yield counts
yield_counts <- list(
  tier_a    = numeric(n_subsets),
  voting_de = numeric(n_subsets)
)
for (q in QUANTIFIERS) yield_counts[[q]] <- numeric(n_subsets)

# Directional split tracking
dir_split_counts <- integer(n_subsets)

for (si in seq_len(n_subsets)) {
  ss <- subsets[[si]]
  tt <- samples_df %>% filter(.data$sample %in% ss, .data$group == CONTRAST_TREAT) %>% pull("sample")
  cc <- samples_df %>% filter(.data$sample %in% ss, .data$group == CONTRAST_CTRL) %>% pull("sample")
  cat(sprintf("\n  ── Subset %d/%d: %s (T) vs %s (C) ──\n",
              si, n_subsets,
              paste(tt, collapse = ", "),
              paste(cc, collapse = ", ")))

  # Run DESeq2 for each quantifier
  data_list <- list()
  for (q in QUANTIFIERS) {
    cat(sprintf("    DESeq2: %s... ", q))
    res <- run_deseq2_subset(count_mats[[q]], ss, q)
    # Filter to universe for consistency
    res <- res %>% filter(.data$gene_id_standard %in% universe)
    data_list[[q]] <- res
    cat(sprintf("%d genes with results\n", nrow(res)))

    # Record single-quantifier DE calls
    de_genes <- res %>%
      filter(!is.na(.data$adj.P.Val), .data$adj.P.Val <= FDR_THRESHOLD,
             abs(.data$logFC) >= LFC_THRESHOLD) %>%
      pull("gene_id_standard")
    genes_by_method[[q]][universe %in% de_genes, si] <- TRUE
    yield_counts[[q]][si] <- length(de_genes)

    # Single-quantifier up-regulated ranking (by padj)
    up_ranks <- res %>%
      filter(!is.na(.data$adj.P.Val), !is.na(.data$logFC),
             .data$logFC > 0) %>%
      arrange(.data$adj.P.Val) %>%
      pull("gene_id_standard")
    # Fill rest of universe at the end
    remaining <- setdiff(universe, up_ranks)
    rank_lists[[q]][[si]] <- c(up_ranks, remaining)
  }

  # Run consensus
  cat("    Consensus: running RRA + voting...\n")
  cons <- run_consensus_subset(data_list, universe)

  # Record Tier A calls
  genes_by_method[["tier_a"]][universe %in% cons$tier_a_genes, si] <- TRUE
  genes_by_method[["voting_de"]][universe %in% cons$voting_de_genes, si] <- TRUE

  rank_lists[["tier_a"]][[si]] <- cons$consensus_rank
  rank_lists[["voting"]][[si]] <- cons$voting_rank

  yield_counts[["tier_a"]][si]    <- cons$n_tier_a
  yield_counts[["voting_de"]][si] <- cons$n_voting_de
  dir_split_counts[si]            <- cons$n_2up_2down

  cat(sprintf("    Tier A: %d | Voting DE: %d | 2-up-2-down splits: %d\n",
              cons$n_tier_a, cons$n_voting_de, cons$n_2up_2down))
}

# ── Step 6: Compute stability metrics ────────────────────────────────────────
cat("\n[6/7] Computing stability metrics...\n")

# Per-method gene-level stability
method_names <- c("tier_a", "voting_de", QUANTIFIERS)
stability_list <- list()
for (m in method_names) {
  stability_list[[m]] <- compute_gene_stability(genes_by_method[[m]])
}

# Build stability table
stability_table <- tibble(
  gene_id = universe,
  stability_consensus    = stability_list[["tier_a"]],
  stability_voting       = stability_list[["voting_de"]],
  stability_featurecounts = stability_list[["featurecounts"]],
  stability_salmon       = stability_list[["salmon"]],
  stability_stringtie    = stability_list[["stringtie"]],
  stability_kallisto     = stability_list[["kallisto"]]
)

# Per-method summary
build_summary_row <- function(name, stab_vec, rank_lst, yield_vec) {
  mean_stab <- mean(stab_vec[stab_vec > 0], na.rm = TRUE)
  if (is.nan(mean_stab) || length(rank_lst) < 2) mean_stab <- NA_real_

  jaccards <- sapply(TOP_N_VALUES, function(n) {
    compute_jaccard_at_n(rank_lst, n)
  })
  names(jaccards) <- paste0("top", TOP_N_VALUES, "_jaccard")

  c(list(method = name, mean_stability = mean_stab,
         total_de_genes = mean(yield_vec, na.rm = TRUE)),
    as.list(jaccards))
}

summary_rows <- list()
for (m in method_names) {
  summary_rows[[m]] <- build_summary_row(
    m, stability_list[[m]], rank_lists[[m]], yield_counts[[m]]
  )
}
summary_table <- bind_rows(lapply(summary_rows, as_tibble))

cat(sprintf("\n  Directional splits (2-up-2-down) per subset:\n"))
cat(sprintf("    Min: %d | Median: %d | Max: %d\n",
            min(dir_split_counts),
            as.integer(median(dir_split_counts)),
            max(dir_split_counts)))

# ── Step 7: Output ───────────────────────────────────────────────────────────
cat("\n[7/7] Writing outputs...\n")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# 7a. Stability TSV
stab_path <- file.path(OUT_DIR, "subsampling_stability.tsv")
write_tsv(stability_table, stab_path)
cat(sprintf("  ✓ %s (%d genes)\n", stab_path, nrow(stability_table)))

# 7b. Summary TSV
sum_path <- file.path(OUT_DIR, "subsampling_summary.tsv")
write_tsv(summary_table, sum_path)
cat(sprintf("  ✓ %s (%d methods)\n", sum_path, nrow(summary_table)))

# 7c. Stability distribution plot
cat("  Generating stability distribution plot...\n")

# Prepare long-format data for plotting
plot_method_labels <- c(
  tier_a = "Consensus Tier A",
  voting_de = "Voting (≥2)",
  featurecounts = "featureCounts",
  salmon = "Salmon",
  stringtie = "StringTie",
  kallisto = "Kallisto"
)

plot_df <- bind_rows(lapply(method_names, function(m) {
  tibble(
    method = factor(plot_method_labels[m], levels = unname(plot_method_labels)),
    stability = stability_list[[m]]
  )
}))

p_dist <- ggplot(plot_df, aes(x = .data$stability, color = .data$method,
                               fill = .data$method)) +
  geom_density(alpha = 0.15, linewidth = 0.8) +
  labs(
    title = "Gene-level DE call stability across 9 subsets",
    subtitle = paste0("Drosophila Wolbachia-infected vs Wolbachia-free | ",
                      length(universe), " common genes | ",
                      "FDR < ", FDR_THRESHOLD, ", |logFC| ≥ ", LFC_THRESHOLD),
    x = "Stability (fraction of subsets where gene is called DE)",
    y = "Density",
    color = "Method", fill = "Method"
  ) +
  theme_bw(base_size = 12) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  coord_cartesian(xlim = c(0.5, 1))

dist_path <- file.path(FIG_DIR, "stability_distribution.png")
ggsave(dist_path, plot = p_dist, width = 9, height = 6, dpi = 150)
cat(sprintf("  ✓ %s\n", dist_path))

# 7d. Stability vs yield curve
cat("  Generating stability vs yield curve...\n")

p_yield <- ggplot(summary_table, aes(x = .data$total_de_genes,
                                      y = .data$mean_stability,
                                      color = .data$method,
                                      label = .data$method)) +
  geom_point(size = 3.5) +
  geom_text(vjust = -1, size = 3.2, show.legend = FALSE) +
  labs(
    title = "Stability vs DEG yield trade-off",
    subtitle = paste0("Drosophila 60d vs 1d | ",
                      length(universe), " common genes | ",
                      "Up-regulated genes only"),
    x = "Mean number of DE genes called (across 9 subsets)",
    y = "Mean gene-level stability",
    color = "Method"
  ) +
  theme_bw(base_size = 12) +
  expand_limits(y = 0)

yield_path <- file.path(FIG_DIR, "stability_yield_curve.png")
ggsave(yield_path, plot = p_yield, width = 9, height = 6, dpi = 150)
cat(sprintf("  ✓ %s\n", yield_path))

# ── Final summary ────────────────────────────────────────────────────────────
cat("\n═══════════════════════════════════════════════════════════════\n")
cat("  BENCHMARK COMPLETE\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat(sprintf("\n  Common universe: %d genes\n", length(universe)))
cat(sprintf("  Subsets: %d (3 choose 2 per group, exhaustive)\n", n_subsets))

cat("\n  Per-method mean yield (± SD across subsets):\n")
for (m in method_names) {
  cat(sprintf("    %-20s  %.1f ± %.1f\n", m,
              mean(yield_counts[[m]], na.rm = TRUE),
              sd(yield_counts[[m]], na.rm = TRUE)))
}

cat(sprintf("\n  Mean stability (genes called in ≥1 subset):\n"))
for (m in method_names) {
  vals <- stability_list[[m]]
  active <- vals[vals > 0]
  if (length(active) > 0) {
    cat(sprintf("    %-25s  %.3f  (n=%d genes with any call)\n",
                m, mean(active), length(active)))
  } else {
    cat(sprintf("    %-25s  N/A (no DE calls)\n", m))
  }
}

cat("\n  Output files:\n")
cat(sprintf("    %s\n", stab_path))
cat(sprintf("    %s\n", sum_path))
cat(sprintf("    %s\n", dist_path))
cat(sprintf("    %s\n", yield_path))
cat("\n")
