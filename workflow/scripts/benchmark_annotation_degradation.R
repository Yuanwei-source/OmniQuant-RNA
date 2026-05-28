#!/usr/bin/env Rscript
# nolint start: object_usage_linter.
# workflow/scripts/benchmark_annotation_degradation.R
# Robustness benchmark: progressively corrupt the GFF annotation file and
# measure how OmniQuant consensus degrades, simulating non-model insect annotations.
#
# Five degradation modes:
#   1. Random drop:  delete X% of gene entries
#   2. Length-biased: delete shortest X% of genes (by genomic span)
#   3. Expression-biased: delete lowest-expressed X% (by mean count)
#   4. Transcript-level: delete X% of transcript isoform/exon entries
#   5. ID corruption: modify X% of gene_id/transcript_id values
#
# Levels: 0% (baseline), 25%, 50%, 75%. 5 random seeds per (mode × level).

# Load packages with box if available, fall back to library()
if (requireNamespace("box", quietly = TRUE)) {
  box::use(
    rd = readr,
    dp = dplyr,
    tb = tibble,
    gg = ggplot2
  )
  `%>%` <- dp$`%>%`
  all_of <- dp$all_of
  any_of <- dp$any_of
  arrange <- dp$arrange
  coalesce <- dp$coalesce
  desc <- dp$desc
  distinct <- dp$distinct
  filter <- dp$filter
  left_join <- dp$left_join
  mutate <- dp$mutate
  select <- dp$select
  transmute <- dp$transmute

  read_csv <- rd$read_csv
  read_tsv <- rd$read_tsv
  write_tsv <- rd$write_tsv

  tibble <- tb$tibble

  aes <- gg$aes
  element_text <- gg$element_text
  element_blank <- gg$element_blank
  geom_line <- gg$geom_line
  geom_point <- gg$geom_point
  geom_ribbon <- gg$geom_ribbon
  ggsave <- gg$ggsave
  ggplot <- gg$ggplot
  labs <- gg$labs
  theme_bw <- gg$theme_bw
  facet_wrap <- gg$facet_wrap
  scale_color_manual <- gg$scale_color_manual
  scale_fill_manual <- gg$scale_fill_manual
} else {
  suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tibble)
    library(ggplot2)
  })
  # dplyr exports via library()
  `%>%` <- dplyr::`%>%`
  all_of <- dplyr::all_of
  any_of <- dplyr::any_of
  arrange <- dplyr::arrange
  coalesce <- dplyr::coalesce
  desc <- dplyr::desc
  distinct <- dplyr::distinct
  filter <- dplyr::filter
  left_join <- dplyr::left_join
  mutate <- dplyr::mutate
  select <- dplyr::select
  transmute <- dplyr::transmute

  read_csv <- readr::read_csv
  read_tsv <- readr::read_tsv
  write_tsv <- readr::write_tsv

  tibble <- tibble::tibble

  aes <- ggplot2::aes
  element_text <- ggplot2::element_text
  geom_line <- ggplot2::geom_line
  geom_point <- ggplot2::geom_point
  geom_ribbon <- ggplot2::geom_ribbon
  ggsave <- ggplot2::ggsave
  ggplot <- ggplot2::ggplot
  labs <- ggplot2::labs
  theme_bw <- ggplot2::theme_bw
  facet_wrap <- ggplot2::facet_wrap
  scale_color_manual <- ggplot2::scale_color_manual
  scale_fill_manual <- ggplot2::scale_fill_manual
}

invisible(utils::globalVariables(c(
  "gene_id_standard", "gene_name", "logFC", "P.Value", "adj.P.Val",
  "included_in_main", "method", "contrast", ".data",
  "mode", "level", "seed", "detectable_recall", "global_recall",
  "precision", "n_detected", "n_reference", "mean_recall", "sd_recall",
  "method_label", "level_pct", "metric", "value"
)))

# ── helpers ────────────────────────────────────────────────────────────────────

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

as_flag <- function(x, default = FALSE) {
  if (is.null(x) || length(x) == 0 || is.na(x[1])) {
    return(default)
  }
  if (is.logical(x)) return(x)
  tolower(as.character(x)) %in% c("true", "1", "yes")
}

clip_probabilities <- function(p, eps) {
  pmin(pmax(p, eps), 1 - eps)
}

coalesce_columns <- function(df, cols) {
  if (length(cols) == 0) return(rep(NA_character_, nrow(df)))
  current <- df[[cols[1]]]
  if (length(cols) == 1) return(current)
  for (col in cols[-1]) {
    current <- dplyr::coalesce(current, df[[col]])
  }
  current
}

# ── consensus functions (replicated from run_consensus_dea.R) ───────────────────

build_rank_list <- function(df, direction) {
  ranked <- df %>%
    filter(!is.na(.data$gene_id_standard), !is.na(.data$P.Value), !is.na(.data$logFC))
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

  result <- tibble(
    gene_id_standard = universe,
    p = rep(1, length(universe))
  )

  if (length(rank_lists) == 0) {
    result$fdr <- p.adjust(result$p, method = "BH")
    return(result)
  }

  rra_out <- RobustRankAggreg::aggregateRanks(
    glist = rank_lists,
    N = length(universe),
    method = "RRA",
    exact = FALSE
  )

  observed <- tibble(
    gene_id_standard = rra_out$Name,
    p_observed = rra_out$Score
  )

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

  tibble(
    gene_id_standard = universe,
    p = cct_p,
    fdr = p.adjust(cct_p, method = "BH")
  )
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

assign_tier <- function(support_n, sign_consistency_n, best_rra_fdr, best_cct_fdr,
                        logfc_cv, tiers, consensus_direction = "up") {
  if (consensus_direction == "mixed" || consensus_direction == "none") {
    return("unclassified")
  }
  if (support_n < 2) return("unclassified")

  if (
    support_n >= tiers$tier_a$min_support &&
    sign_consistency_n >= tiers$tier_a$min_sign_consistency &&
    sign_consistency_n == support_n &&
    best_rra_fdr <= tiers$tier_a$max_rra_fdr &&
    best_cct_fdr <= tiers$tier_a$max_cct_fdr &&
    (is.na(logfc_cv) || logfc_cv <= tiers$tier_a$max_logfc_cv)
  ) return("Tier_A")

  if (
    support_n >= tiers$tier_b$min_support &&
    sign_consistency_n >= tiers$tier_b$min_sign_consistency &&
    best_rra_fdr <= tiers$tier_b$max_rra_fdr &&
    best_cct_fdr <= tiers$tier_b$max_cct_fdr &&
    (is.na(logfc_cv) || logfc_cv <= tiers$tier_b$max_logfc_cv)
  ) return("Tier_B")

  if (
    support_n >= tiers$tier_c$min_support &&
    sign_consistency_n >= tiers$tier_c$min_sign_consistency &&
    best_rra_fdr <= tiers$tier_c$max_rra_fdr &&
    best_cct_fdr <= tiers$tier_c$max_cct_fdr &&
    (is.na(logfc_cv) || logfc_cv <= tiers$tier_c$max_logfc_cv)
  ) return("Tier_C")

  "unclassified"
}

compute_logfc_cv <- function(logfc_matrix, consensus_direction) {
  res <- numeric(nrow(logfc_matrix))
  for (i in seq_len(nrow(logfc_matrix))) {
    values <- logfc_matrix[i, ]
    direction <- consensus_direction[i]
    if (direction == "mixed" || direction == "none") {
      res[i] <- NA_real_
      next
    }
    values <- values[!is.na(values)]
    if (length(values) < 2) { res[i] <- NA_real_; next }
    if (direction == "up")      values <- values[values > 0]
    else if (direction == "down") values <- abs(values[values < 0])
    else values <- abs(values)
    if (length(values) < 2) { res[i] <- NA_real_; next }
    mean_value <- mean(values)
    if (is.na(mean_value) || mean_value == 0) { res[i] <- NA_real_; next }
    res[i] <- sd(values) / mean_value
  }
  res
}

# ── run full consensus pipeline (adapted from run_consensus_dea.R) ──────────────

run_consensus <- function(data_list, universe, quantifier_names,
                          fdr_threshold, lfc_threshold, p_clip, tiers) {
  consensus_df <- tibble(gene_id_standard = universe)

  # Build joined data frame
  for (qn in quantifier_names) {
    df <- data_list[[qn]] %>%
      transmute(
        gene_id_standard = .data$gene_id_standard,
        !!paste0("logFC__", qn) := .data$logFC,
        !!paste0("P.Value__", qn) := .data$P.Value,
        !!paste0("adj.P.Val__", qn) := .data$adj.P.Val
      )
    consensus_df <- consensus_df %>% left_join(df, by = "gene_id_standard")
  }

  logfc_cols <- paste0("logFC__", quantifier_names)
  p_cols     <- paste0("P.Value__", quantifier_names)
  adj_cols   <- paste0("adj.P.Val__", quantifier_names)

  logfc_matrix <- as.matrix(consensus_df[, logfc_cols, drop = FALSE])
  p_matrix     <- as.matrix(consensus_df[, p_cols, drop = FALSE])
  adj_matrix   <- as.matrix(consensus_df[, adj_cols, drop = FALSE])
  storage.mode(logfc_matrix) <- "numeric"
  storage.mode(p_matrix)     <- "numeric"
  storage.mode(adj_matrix)   <- "numeric"

  sig_matrix <- !is.na(logfc_matrix) & !is.na(adj_matrix) &
    (adj_matrix <= fdr_threshold) & (abs(logfc_matrix) >= lfc_threshold)

  support_n       <- rowSums(sig_matrix, na.rm = TRUE)
  up_support_n    <- rowSums(sig_matrix & logfc_matrix > 0, na.rm = TRUE)
  down_support_n  <- rowSums(sig_matrix & logfc_matrix < 0, na.rm = TRUE)
  sign_consistency_n <- pmax(up_support_n, down_support_n)

  rra_up   <- compute_rra_scores(data_list, universe, direction = "up")
  rra_down <- compute_rra_scores(data_list, universe, direction = "down")
  colnames(rra_up)   <- c("gene_id_standard", "rra_up_p", "rra_up_fdr")
  colnames(rra_down) <- c("gene_id_standard", "rra_down_p", "rra_down_fdr")

  cct_up   <- compute_cct_scores(logfc_matrix, p_matrix, universe, "up", p_clip)
  cct_down <- compute_cct_scores(logfc_matrix, p_matrix, universe, "down", p_clip)
  colnames(cct_up)   <- c("gene_id_standard", "cct_up_p", "cct_up_fdr")
  colnames(cct_down) <- c("gene_id_standard", "cct_down_p", "cct_down_fdr")

  consensus_df <- consensus_df %>%
    left_join(rra_up, by = "gene_id_standard") %>%
    left_join(rra_down, by = "gene_id_standard") %>%
    left_join(cct_up, by = "gene_id_standard") %>%
    left_join(cct_down, by = "gene_id_standard")

  consensus_direction <- vapply(
    seq_len(nrow(consensus_df)),
    function(i) choose_direction(up_support_n[i], down_support_n[i]),
    character(1)
  )

  logfc_cv <- compute_logfc_cv(logfc_matrix, consensus_direction)

  best_rra_p <- ifelse(consensus_direction == "up",
    consensus_df$rra_up_p,
    ifelse(consensus_direction == "down", consensus_df$rra_down_p, NA_real_))
  best_rra_fdr <- ifelse(consensus_direction == "up",
    consensus_df$rra_up_fdr,
    ifelse(consensus_direction == "down", consensus_df$rra_down_fdr, NA_real_))
  best_cct_fdr <- ifelse(consensus_direction == "up",
    consensus_df$cct_up_fdr,
    ifelse(consensus_direction == "down", consensus_df$cct_down_fdr, NA_real_))

  tier <- vapply(
    seq_len(nrow(consensus_df)),
    function(i) {
      assign_tier(
        support_n         = support_n[i],
        sign_consistency_n = sign_consistency_n[i],
        best_rra_fdr       = best_rra_fdr[i],
        best_cct_fdr       = best_cct_fdr[i],
        logfc_cv           = logfc_cv[i],
        tiers              = tiers,
        consensus_direction = consensus_direction[i]
      )
    },
    character(1)
  )

  # Per-quantifier significance
  per_quant_sig <- as.data.frame(sig_matrix)
  colnames(per_quant_sig) <- quantifier_names

  list(
    df        = consensus_df,
    tier      = tier,
    sig_matrix = sig_matrix,
    support_n = support_n,
    per_quant_sig = per_quant_sig
  )
}

# ── GFF parsing ─────────────────────────────────────────────────────────────────

parse_gff <- function(gff_path) {
  cat("Parsing GFF file:", gff_path, "\n")
  lines <- readLines(gff_path)
  data_lines <- lines[!grepl("^#", lines) & nchar(trimws(lines)) > 0]
  cat("  Non-comment lines:", length(data_lines), "\n")

  parts <- strsplit(data_lines, "\t")
  valid <- lengths(parts) >= 9
  parts <- parts[valid]

  seqid     <- vapply(parts, `[`, character(1), 1)
  source    <- vapply(parts, `[`, character(1), 2)
  feature   <- vapply(parts, `[`, character(1), 3)
  start     <- as.integer(vapply(parts, `[`, character(1), 4))
  end       <- as.integer(vapply(parts, `[`, character(1), 5))
  strand    <- vapply(parts, `[`, character(1), 7)
  attributes_str <- vapply(parts, `[`, character(1), 9)

  # Extract key attributes
  extract_attr <- function(attr_str, key) {
    pattern <- paste0(key, "=([^;]+)")
    val <- regmatches(attr_str, regexec(pattern, attr_str))
    vapply(val, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1))
  }

  gene_id       <- extract_attr(attributes_str, "gene_id")
  transcript_id <- extract_attr(attributes_str, "transcript_id")
  parent        <- extract_attr(attributes_str, "Parent")
  biotype       <- extract_attr(attributes_str, "biotype")
  id_attr       <- extract_attr(attributes_str, "ID")  # Full ID with prefix (gene:FBgnXXXX)

  # For gene entries, use the ID attribute (includes "gene:" prefix) to match
  # the namespace used by featureCounts/consensus (gene:FBgnXXXX format).
  # For non-gene entries, fall back to gene_id attribute.
  # Standardize: prepend "gene:" to gene_id if not already prefixed
  gene_id_std <- ifelse(
    !is.na(id_attr) & grepl("^gene:", id_attr),
    id_attr,
    ifelse(!is.na(gene_id) & !grepl("^gene:", gene_id),
           paste0("gene:", gene_id),
           gene_id)
  )

  # Create data frame using base R (avoid tibble column name issues)
  df <- data.frame(
    seqid     = seqid,
    source    = source,
    feature   = feature,
    start     = start,
    end       = end,
    strand    = strand,
    gene_id   = gene_id_std,   # standardized: uses "gene:FBgnXXXX" format
    gene_id_raw = gene_id,     # raw attribute value (FBgnXXXX)
    transcript_id = transcript_id,
    parent    = parent,
    biotype   = biotype,
    stringsAsFactors = FALSE
  )

  cat("  Parsed features:", nrow(df), "\n")
  cat("  Feature types:", paste(sort(unique(df$feature)), collapse = ", "), "\n")
  df
}

# ── degradation mode functions ─────────────────────────────────────────────────

degrade_random_drop <- function(gff_df, gene_genes, level, seed) {
  # level is fraction 0-1
  set.seed(seed)
  all_genes <- unique(gene_genes$gene_id)
  all_genes <- all_genes[!is.na(all_genes) & all_genes != ""]
  cat(sprintf("  [random_drop] total genes: %d\n", length(all_genes)))

  if (level == 0) {
    return(list(dropped_genes = character(0), retained_genes = all_genes))
  }

  n_drop <- round(length(all_genes) * level)
  dropped <- sample(all_genes, n_drop)
  retained <- setdiff(all_genes, dropped)
  cat(sprintf("  [random_drop] level=%.2f: dropped %d, retained %d\n",
              level, length(dropped), length(retained)))
  list(dropped_genes = dropped, retained_genes = retained)
}

degrade_length_biased <- function(gff_df, gene_genes, level, seed) {
  set.seed(seed)
  # Compute genomic span for each gene
  gene_spans <- gene_genes %>%
    filter(!is.na(.data$gene_id), .data$gene_id != "") %>%
    mutate(span = .data$end - .data$start + 1)

  if (nrow(gene_spans) == 0) {
    return(list(dropped_genes = character(0), retained_genes = character(0)))
  }

  sorted_genes <- gene_spans[order(gene_spans$span), ]
  all_genes <- sorted_genes$gene_id

  if (level == 0) {
    return(list(dropped_genes = character(0), retained_genes = all_genes))
  }

  n_drop <- round(length(all_genes) * level)
  dropped <- all_genes[seq_len(n_drop)]
  retained <- setdiff(all_genes, dropped)
  cat(sprintf("  [length_biased] level=%.2f: dropped %d (shortest), retained %d\n",
              level, length(dropped), length(retained)))
  list(dropped_genes = dropped, retained_genes = retained)
}

degrade_expression_biased <- function(gff_df, gene_genes, level, seed, expr_vec) {
  set.seed(seed)
  # expr_vec: named numeric vector of mean expression (names = gene_ids)
  all_genes <- unique(gene_genes$gene_id)
  all_genes <- all_genes[!is.na(all_genes) & all_genes != ""]

  if (level == 0) {
    return(list(dropped_genes = character(0), retained_genes = all_genes))
  }

  # Get expression for genes present in GFF
  common_genes <- intersect(all_genes, names(expr_vec))
  expr_sub <- expr_vec[common_genes]
  sorted_genes <- names(sort(expr_sub))

  n_drop <- round(length(all_genes) * level)
  n_drop <- min(n_drop, length(sorted_genes))

  dropped <- sorted_genes[seq_len(n_drop)]
  # Also drop genes not in expr_vec (no expression data → can't assess)
  extra_dropped <- setdiff(all_genes, names(expr_vec))
  dropped <- unique(c(dropped, extra_dropped))
  retained <- setdiff(all_genes, dropped)

  cat(sprintf("  [expr_biased] level=%.2f: dropped %d (lowest expr), retained %d\n",
              level, length(dropped), length(retained)))
  list(dropped_genes = dropped, retained_genes = retained)
}

degrade_transcript_level <- function(gff_df, level, seed) {
  set.seed(seed)
  # Target: mRNA, transcript, ncRNA, pre_miRNA, miRNA, etc. — transcript-like features
  transcript_types <- c("mRNA", "transcript", "ncRNA", "pre_miRNA", "miRNA",
                        "tRNA", "snoRNA", "snRNA", "rRNA",
                        "pseudogenic_transcript", "transposable_element")
  transcript_rows <- which(gff_df$feature %in% transcript_types)

  cat(sprintf("  [transcript_level] total transcript entries: %d\n", length(transcript_rows)))

  if (level == 0 || length(transcript_rows) == 0) {
    return(list(dropped_rows = integer(0), affected_genes = character(0),
                dropped_genes = character(0)))
  }

  n_drop <- round(length(transcript_rows) * level)
  drop_idx <- sample(transcript_rows, n_drop)

  # For each dropped transcript, also drop its child exons/CDS/UTR entries
  # Children reference parent via "Parent=transcript:XXX"
  dropped_tx_ids <- unique(gff_df$transcript_id[drop_idx])
  dropped_tx_ids <- dropped_tx_ids[!is.na(dropped_tx_ids) & dropped_tx_ids != ""]

  child_rows <- which(gff_df$parent %in% paste0("transcript:", dropped_tx_ids))
  all_dropped_rows <- unique(c(drop_idx, child_rows))

  # Determine which genes lost ALL of their transcripts
  # For each gene, count remaining transcript entries
  remaining_gff <- gff_df[-all_dropped_rows, ]
  remaining_transcripts <- remaining_gff[
    remaining_gff$feature %in% transcript_types, ]
  remaining_tx_per_gene <- table(remaining_transcripts$gene_id)

  # Genes that now have 0 transcript entries
  original_genes <- unique(gff_df$gene_id[gff_df$feature %in% transcript_types])
  original_genes <- original_genes[!is.na(original_genes) & original_genes != ""]
  dropped_genes <- setdiff(
    original_genes,
    names(remaining_tx_per_gene)[remaining_tx_per_gene > 0]
  )

  cat(sprintf("  [transcript_level] level=%.2f: dropped %d transcript rows, %d genes lost all isoforms\n",
              level, length(all_dropped_rows), length(dropped_genes)))

  list(
    dropped_rows   = all_dropped_rows,
    affected_genes = character(0),
    dropped_genes  = dropped_genes
  )
}

degrade_id_corruption <- function(gff_df, gene_genes, level, seed) {
  set.seed(seed)
  all_genes <- unique(gene_genes$gene_id)
  all_genes <- all_genes[!is.na(all_genes) & all_genes != ""]

  if (level == 0) {
    return(list(corrupted_genes = character(0),
                original_ids     = character(0),
                corrupted_ids    = character(0)))
  }

  n_corrupt <- round(length(all_genes) * level)
  corrupt_targets <- sample(all_genes, n_corrupt)

  # Create corrupted IDs: append random hex suffix
  rand_suffix <- sprintf("%04x", sample(0:65535, length(corrupt_targets), replace = TRUE))
  corrupted <- paste0(corrupt_targets, "_CORRUPT_", rand_suffix)

  cat(sprintf("  [id_corruption] level=%.2f: corrupted %d gene IDs\n",
              level, length(corrupt_targets)))

  list(
    corrupted_genes = corrupt_targets,
    original_ids    = corrupt_targets,
    corrupted_ids   = corrupted
  )
}

# ── apply degradation to quantifier DEA results ─────────────────────────────────

apply_gene_drop_to_dea <- function(dea_df, dropped_genes, quantifier) {
  # For annotation-dependent quantifiers (featurecounts, stringtie):
  #   genes in dropped_genes → set logFC=NA, P.Value=1, adj.P.Val=1
  # For annotation-independent (salmon, kallisto): no change
  annotation_dependent <- c("featurecounts", "stringtie")

  if (!(quantifier %in% annotation_dependent)) {
    return(dea_df)
  }

  if (length(dropped_genes) == 0) return(dea_df)

  # Use base R column assignment (NOT dplyr::rename)
  colnames_dea <- colnames(dea_df)
  id_col <- if ("gene_id_standard" %in% colnames_dea) "gene_id_standard" else "gene_id"

  mask <- dea_df[[id_col]] %in% dropped_genes
  if (sum(mask) > 0) {
    dea_df[mask, "logFC"]    <- NA_real_
    dea_df[mask, "P.Value"]  <- 1
    dea_df[mask, "adj.P.Val"] <- 1
  }

  dea_df
}

apply_id_corruption_to_dea <- function(dea_df, corrupted_genes, corrupted_ids,
                                        quantifier) {
  # For ALL quantifiers: corrupted gene_ids → change to corrupted version
  # This simulates namespace mapping failure
  if (length(corrupted_genes) == 0) return(dea_df)

  id_col <- if ("gene_id_standard" %in% colnames(dea_df)) "gene_id_standard" else "gene_id"
  mapping <- setNames(corrupted_ids, corrupted_genes)

  mask <- dea_df[[id_col]] %in% corrupted_genes
  if (sum(mask) > 0) {
    ids_to_replace <- dea_df[[id_col]][mask]
    dea_df[mask, id_col] <- mapping[ids_to_replace]
  }

  dea_df
}

# ── compute benchmark metrics ───────────────────────────────────────────────────

compute_metrics <- function(result_list, reference_genes, retained_genes) {
  # result_list: output of run_consensus()
  # reference_genes: full-annotation Tier A genes (gold standard)
  # retained_genes: genes still present in degraded annotation universe

  tier_labels <- result_list$tier
  gene_ids    <- result_list$df$gene_id_standard

  # Consensus detection: Tier_A from RRA-driven consensus
  consensus_detected <- gene_ids[tier_labels == "Tier_A"]

  # Per-quantifier detection
  quant_names <- colnames(result_list$per_quant_sig)
  per_quant_detected <- lapply(quant_names, function(qn) {
    gene_ids[result_list$per_quant_sig[[qn]]]
  })
  names(per_quant_detected) <- quant_names

  # Detectable recall: only among genes still in degraded universe
  detectable_universe <- intersect(reference_genes, retained_genes)
  n_ref_detectable   <- length(detectable_universe)

  if (n_ref_detectable > 0) {
    consensus_detectable_recall <- length(intersect(consensus_detected, detectable_universe)) /
                                   n_ref_detectable
    per_q_detectable <- vapply(per_quant_detected, function(det) {
      length(intersect(det, detectable_universe)) / n_ref_detectable
    }, numeric(1))
  } else {
    consensus_detectable_recall <- NA_real_
    per_q_detectable <- rep(NA_real_, length(quant_names))
  }
  names(per_q_detectable) <- paste0("detectable_recall__", quant_names)

  # Global recovery recall: among ALL reference genes
  n_ref_total <- length(reference_genes)

  if (n_ref_total > 0) {
    consensus_global_recall <- length(intersect(consensus_detected, reference_genes)) /
                               n_ref_total
    per_q_global <- vapply(per_quant_detected, function(det) {
      length(intersect(det, reference_genes)) / n_ref_total
    }, numeric(1))
  } else {
    consensus_global_recall <- NA_real_
    per_q_global <- rep(NA_real_, length(quant_names))
  }
  names(per_q_global) <- paste0("global_recall__", quant_names)

  # Precision: fraction of detected DEGs that are in reference set
  n_consensus_detected <- length(consensus_detected)
  if (n_consensus_detected > 0) {
    consensus_precision <- length(intersect(consensus_detected, reference_genes)) /
                           n_consensus_detected
  } else {
    consensus_precision <- NA_real_
  }

  per_q_precision <- vapply(per_quant_detected, function(det) {
    n_det <- length(det)
    if (n_det > 0) length(intersect(det, reference_genes)) / n_det else NA_real_
  }, numeric(1))
  names(per_q_precision) <- paste0("precision__", quant_names)

  c(
    consensus_detectable_recall = consensus_detectable_recall,
    consensus_global_recall     = consensus_global_recall,
    consensus_precision         = consensus_precision,
    n_detected                  = n_consensus_detected,
    n_reference                 = n_ref_total,
    per_q_detectable,
    per_q_global,
    per_q_precision
  )
}

# ── main ────────────────────────────────────────────────────────────────────────

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  # Default paths (can be overridden)
  gff_path  <- "data/reference/Drosophila_melanogaster.BDGP6.54.115.chr.gff3"
  fc_dea    <- "results/06.differential_expression/featurecounts/deseq2.60d_vs_1d.csv"
  st_dea    <- "results/06.differential_expression/stringtie/deseq2.60d_vs_1d.csv"
  sa_dea    <- "results/06.differential_expression/salmon/deseq2.60d_vs_1d.csv"
  ka_dea    <- "results/06.differential_expression/kallisto/deseq2.60d_vs_1d.csv"
  fc_counts <- "results/05.quantification/matrices/featurecounts/featurecounts_gene_counts_matrix.tsv"
  output_tsv <- "results/benchmark/annotation_degradation.tsv"
  output_fig <- "results/benchmark/figures/degradation_curve.png"
  n_seeds    <- 5L
  contrast   <- "Wolbachia_infected_vs_Wolbachia_free"
  levels_vec <- c(0.0, 0.25, 0.50, 0.75)

  # Help message
  if (length(args) == 0 && !interactive()) {
    cat("Usage: Rscript workflow/scripts/benchmark_annotation_degradation.R [OPTIONS]\n\n")
    cat("Options:\n")
    cat("  --gff PATH               GFF annotation file\n")
    cat("  --featurecounts-dea PATH featureCounts DEA CSV\n")
    cat("  --stringtie-dea PATH     StringTie DEA CSV\n")
    cat("  --salmon-dea PATH        Salmon DEA CSV\n")
    cat("  --kallisto-dea PATH      Kallisto DEA CSV\n")
    cat("  --featurecounts-counts PATH featureCounts count matrix\n")
    cat("  --output-tsv PATH        Output per-seed metrics TSV\n")
    cat("  --output-figure PATH     Output degradation curve PNG\n")
    cat("  --seeds N                Random seeds per (mode x level) [5]\n")
    cat("  --contrast NAME          Contrast name [60d_vs_1d]\n")
    cat("  --levels 0,0.25,0.5,0.75  Degradation levels as fractions or % (comma-sep) [0,0.25,0.5,0.75]\n\n")
    cat("Requires: R packages dplyr, readr, tibble, ggplot2, RobustRankAggreg\n")
    cat("          (install via: conda env create -f envs/dea.yaml)\n")
    return(invisible(NULL))
  }

  # Parse args
  i <- 1
  while (i <= length(args)) {
    switch(args[i],
      "--gff"               = { gff_path  <- args[i + 1]; i <- i + 2 },
      "--featurecounts-dea" = { fc_dea    <- args[i + 1]; i <- i + 2 },
      "--stringtie-dea"     = { st_dea    <- args[i + 1]; i <- i + 2 },
      "--salmon-dea"        = { sa_dea    <- args[i + 1]; i <- i + 2 },
      "--kallisto-dea"      = { ka_dea    <- args[i + 1]; i <- i + 2 },
      "--featurecounts-counts" = { fc_counts <- args[i + 1]; i <- i + 2 },
      "--output-tsv"        = { output_tsv <- args[i + 1]; i <- i + 2 },
      "--output-figure"     = { output_fig <- args[i + 1]; i <- i + 2 },
      "--seeds"             = { n_seeds   <- as.integer(args[i + 1]); i <- i + 2 },
      "--contrast"          = { contrast  <- args[i + 1]; i <- i + 2 },
      "--levels"            = { levels_vec <- as.numeric(strsplit(args[i + 1], ",")[[1]]); i <- i + 2 },
      { cat("Unknown argument:", args[i], "\n"); quit(status = 1) }
    )
  }

  # Normalize levels: if any value > 1, treat as percentages and divide by 100
  if (any(levels_vec > 1)) {
    levels_vec <- levels_vec / 100
  }
  levels_vec <- sort(unique(levels_vec))

  cat("\n========== Annotation Degradation Benchmark ==========\n")
  cat("GFF:", gff_path, "\n")
  cat("Contrast:", contrast, "\n")
  cat("Levels:", paste(levels_vec, collapse = ", "), "\n")
  cat("Seeds per (mode × level):", n_seeds, "\n\n")

  # ── tier configuration (CCT fdr = 1.0, RRA-driven consensus only) ──
  tiers <- list(
    tier_a = list(
      min_support = 4L, min_sign_consistency = 4L,
      max_rra_fdr = 0.05, max_cct_fdr = 1.0, max_logfc_cv = 1.0
    ),
    tier_b = list(
      min_support = 3L, min_sign_consistency = 3L,
      max_rra_fdr = 0.10, max_cct_fdr = 1.0, max_logfc_cv = 1.25
    ),
    tier_c = list(
      min_support = 2L, min_sign_consistency = 2L,
      max_rra_fdr = 0.25, max_cct_fdr = 1.0, max_logfc_cv = 1.50
    )
  )
  fdr_threshold <- 0.05
  lfc_threshold <- 1.0
  p_clip <- 1e-16
  quantifier_names <- c("featurecounts", "stringtie", "salmon", "kallisto")

  # ── load DEA results ─────────────────────────────────────
  cat("── Loading DEA results ──\n")

  load_dea <- function(path, quantifier, contrast_name) {
    df <- read_csv(path, show_col_types = FALSE) %>%
      filter(.data$contrast == contrast_name,
             as_flag(.data$included_in_main, default = TRUE)) %>%
      mutate(
        quantifier = quantifier,
        gene_id_standard = coalesce(.data$gene_id_standard, .data$gene_id)
      ) %>%
      filter(!is.na(.data$gene_id_standard), .data$gene_id_standard != "") %>%
      distinct(.data$gene_id_standard, .keep_all = TRUE) %>%
      select("gene_id_standard", "gene_name", "logFC", "P.Value", "adj.P.Val",
             "quantifier")
    df
  }

  fc_res <- load_dea(fc_dea, "featurecounts", contrast)
  st_res <- load_dea(st_dea, "stringtie", contrast)
  sa_res <- load_dea(sa_dea, "salmon", contrast)
  ka_res <- load_dea(ka_dea, "kallisto", contrast)

  cat(sprintf("  featurecounts: %d genes\n", nrow(fc_res)))
  cat(sprintf("  stringtie:     %d genes\n", nrow(st_res)))
  cat(sprintf("  salmon:        %d genes\n", nrow(sa_res)))
  cat(sprintf("  kallisto:      %d genes\n", nrow(ka_res)))

  # ── build universe (union of all quantifier genes) ────────
  all_genes <- unique(c(
    fc_res$gene_id_standard, st_res$gene_id_standard,
    sa_res$gene_id_standard, ka_res$gene_id_standard
  ))
  universe <- sort(all_genes)
  cat(sprintf("\n  Universe: %d genes\n", length(universe)))

  data_list <- list(
    featurecounts = fc_res,
    stringtie     = st_res,
    salmon        = sa_res,
    kallisto      = ka_res
  )

  # ── baseline consensus (full annotation, 0% degradation) ─
  cat("\n── Running baseline consensus ──\n")
  baseline <- run_consensus(
    data_list, universe, quantifier_names,
    fdr_threshold, lfc_threshold, p_clip, tiers
  )

  reference_genes <- universe[baseline$tier == "Tier_A"]
  cat(sprintf("  Reference set (full-annotation Tier A): %d genes\n",
              length(reference_genes)))

  # ── parse GFF ──────────────────────────────────────────────
  cat("\n── Parsing GFF ──\n")
  gff_df <- parse_gff(gff_path)

  # Extract gene entries (feature == "gene")
  gene_genes <- gff_df[gff_df$feature == "gene", ]
  cat(sprintf("  Gene entries: %d\n", nrow(gene_genes)))

  # ── compute expression vector for mode 3 ──────────────────
  cat("\n── Loading featureCounts counts for expression-biased degradation ──\n")
  fc_counts_df <- read_tsv(fc_counts, comment = "#", show_col_types = FALSE)
  # Standardize first column name
  colnames(fc_counts_df)[1] <- "gene_id"
  # Compute mean expression per gene (exclude annotation columns if present)
  count_cols <- setdiff(colnames(fc_counts_df),
                        c("gene_id", "Chr", "Start", "End", "Strand", "Length"))
  count_mat <- as.matrix(fc_counts_df[, count_cols, drop = FALSE])
  storage.mode(count_mat) <- "numeric"
  expr_vec <- rowMeans(count_mat, na.rm = TRUE)
  names(expr_vec) <- fc_counts_df$gene_id
  cat(sprintf("  Expression vector: %d genes\n", length(expr_vec)))

  # ── define degradation modes ──────────────────────────────
  modes <- list(
    list(name = "random_drop",       fn = degrade_random_drop,
         desc = "Random gene drop"),
    list(name = "length_biased",     fn = degrade_length_biased,
         desc = "Length-biased gene drop"),
    list(name = "expression_biased", fn = degrade_expression_biased,
         desc = "Expression-biased gene drop"),
    list(name = "transcript_level",  fn = NULL,
         desc = "Transcript-level entry drop")
  )

  # ── benchmark loop ────────────────────────────────────────
  cat("\n── Running benchmark grid ──\n")
  results_list <- list()
  total_runs <- length(modes) * length(levels_vec) * n_seeds
  run_count <- 0

  for (m in seq_along(modes)) {
    mode_info <- modes[[m]]
    cat(sprintf("\n--- Mode: %s (%s) ---\n", mode_info$name, mode_info$desc))

    for (level in levels_vec) {
      cat(sprintf("  Level: %.0f%%\n", level * 100))

      for (seed in seq_len(n_seeds)) {
        run_count <- run_count + 1
        if (run_count %% 10 == 0) {
          cat(sprintf("    ... %d/%d runs completed\n", run_count, total_runs))
        }

        actual_seed <- seed * 100 + m  # unique per mode×level×seed

        # Apply degradation
        if (mode_info$name == "transcript_level") {
          deg <- degrade_transcript_level(gff_df, level, actual_seed)
          dropped_genes <- deg$dropped_genes
          affected_genes <- deg$affected_genes
          # For transcript-level: genes that lost all isoforms → dropped
          # Use same logic as gene-level drop for annotation-dependent quantifiers
        } else if (mode_info$name == "id_corruption") {
          deg <- degrade_id_corruption(gff_df, gene_genes, level, actual_seed)
          # ID corruption: ALL quantifiers affected (namespace mismatch)
        } else {
          if (mode_info$name == "expression_biased") {
            deg <- mode_info$fn(gff_df, gene_genes, level, actual_seed, expr_vec)
          } else {
            deg <- mode_info$fn(gff_df, gene_genes, level, actual_seed)
          }
        }

        # Apply to quantifier data
        if (mode_info$name == "id_corruption") {
          # ID corruption: change IDs in all quantifiers
          mod_fc <- apply_id_corruption_to_dea(
            fc_res, deg$corrupted_genes, deg$corrupted_ids, "featurecounts")
          mod_st <- apply_id_corruption_to_dea(
            st_res, deg$corrupted_genes, deg$corrupted_ids, "stringtie")
          mod_sa <- apply_id_corruption_to_dea(
            sa_res, deg$corrupted_genes, deg$corrupted_ids, "salmon")
          mod_ka <- apply_id_corruption_to_dea(
            ka_res, deg$corrupted_genes, deg$corrupted_ids, "kallisto")
          # Retained genes = genes NOT corrupted
          corrupted_set <- deg$corrupted_genes
          all_gff_genes <- unique(gene_genes$gene_id)
          all_gff_genes <- all_gff_genes[!is.na(all_gff_genes) & all_gff_genes != ""]
          retained_genes <- setdiff(all_gff_genes, corrupted_set)
          dropped_genes <- corrupted_set
        } else if (mode_info$name == "transcript_level") {
          mod_fc <- apply_gene_drop_to_dea(fc_res, deg$dropped_genes, "featurecounts")
          mod_st <- apply_gene_drop_to_dea(st_res, deg$dropped_genes, "stringtie")
          mod_sa <- sa_res  # annotation-independent
          mod_ka <- ka_res  # annotation-independent
          retained_genes <- setdiff(
            unique(gene_genes$gene_id[!is.na(gene_genes$gene_id) & gene_genes$gene_id != ""]),
            deg$dropped_genes
          )
          dropped_genes <- deg$dropped_genes
        } else {
          # Gene-level drops: affects featurecounts and stringtie only
          mod_fc <- apply_gene_drop_to_dea(fc_res, deg$dropped_genes, "featurecounts")
          mod_st <- apply_gene_drop_to_dea(st_res, deg$dropped_genes, "stringtie")
          mod_sa <- sa_res  # unchanged
          mod_ka <- ka_res  # unchanged
          retained_genes <- deg$retained_genes
          dropped_genes <- deg$dropped_genes
        }

        # Rebuild universe with modified data
        mod_universe <- sort(unique(c(
          mod_fc$gene_id_standard, mod_st$gene_id_standard,
          mod_sa$gene_id_standard, mod_ka$gene_id_standard
        )))

        mod_data <- list(
          featurecounts = mod_fc,
          stringtie     = mod_st,
          salmon        = mod_sa,
          kallisto      = mod_ka
        )

        # Run consensus
        result <- tryCatch(
          run_consensus(mod_data, mod_universe, quantifier_names,
                        fdr_threshold, lfc_threshold, p_clip, tiers),
          error = function(e) {
            cat(sprintf("      ERROR in consensus: %s\n", conditionMessage(e)))
            NULL
          }
        )

        if (is.null(result)) next

        # Compute metrics
        metrics <- compute_metrics(result, reference_genes, retained_genes)

        # Consensus row (once per run)
        row_consensus <- data.frame(
          mode    = mode_info$name,
          level   = level,
          method  = "consensus",
          detectable_recall = as.numeric(metrics["consensus_detectable_recall"]),
          global_recall     = as.numeric(metrics["consensus_global_recall"]),
          precision         = as.numeric(metrics["consensus_precision"]),
          n_detected        = as.integer(metrics["n_detected"]),
          n_reference       = as.integer(metrics["n_reference"]),
          seed    = seed,
          stringsAsFactors = FALSE
        )
        results_list[[length(results_list) + 1]] <- row_consensus

        # Per-quantifier rows
        for (qn in quantifier_names) {
          dr_key <- paste0("detectable_recall__", qn)
          gr_key <- paste0("global_recall__", qn)
          pr_key <- paste0("precision__", qn)

          row_q <- data.frame(
            mode    = mode_info$name,
            level   = level,
            method  = qn,
            detectable_recall = as.numeric(metrics[dr_key]),
            global_recall     = as.numeric(metrics[gr_key]),
            precision         = as.numeric(metrics[pr_key]),
            n_detected        = NA_integer_,
            n_reference       = as.integer(metrics["n_reference"]),
            seed    = seed,
            stringsAsFactors = FALSE
          )
          results_list[[length(results_list) + 1]] <- row_q
        }
      }
    }
  }

  # ── aggregate results ─────────────────────────────────────
  cat("\n── Aggregating results ──\n")
  all_results <- do.call(rbind, results_list)
  rownames(all_results) <- NULL

  # Compute mean ± SD across seeds
  aggr <- all_results %>%
    dplyr::group_by(mode, level, method) %>%
    dplyr::summarise(
      mean_detectable_recall = mean(detectable_recall, na.rm = TRUE),
      sd_detectable_recall   = sd(detectable_recall, na.rm = TRUE),
      mean_global_recall     = mean(global_recall, na.rm = TRUE),
      sd_global_recall       = sd(global_recall, na.rm = TRUE),
      mean_precision         = mean(precision, na.rm = TRUE),
      sd_precision           = sd(precision, na.rm = TRUE),
      mean_n_detected        = mean(n_detected, na.rm = TRUE),
      n_seeds                = dplyr::n(),
      .groups = "drop"
    )

  # Write raw per-seed TSV
  dir.create(dirname(output_tsv), recursive = TRUE, showWarnings = FALSE)
  write_tsv(all_results, output_tsv)
  cat(sprintf("  Wrote %d rows to %s\n", nrow(all_results), output_tsv))

  # Write aggregated summary
  aggr_path <- file.path(dirname(output_tsv), "annotation_degradation_summary.tsv")
  write_tsv(aggr, aggr_path)
  cat(sprintf("  Wrote summary to %s\n", aggr_path))

  # ── generate degradation curve figure ─────────────────────
  cat("\n── Generating figures ──\n")
  dir.create(dirname(output_fig), recursive = TRUE, showWarnings = FALSE)

  # Prepare plot data
  plot_df <- aggr %>%
    mutate(
      level_pct = level * 100,
      method_label = method
    )

  # Color palette
  method_colors <- c(
    "consensus"      = "#000000",
    "featurecounts"  = "#E41A1C",
    "stringtie"      = "#377EB8",
    "salmon"         = "#4DAF4A",
    "kallisto"       = "#984EA3"
  )

  mode_labels <- c(
    "random_drop"       = "Random drop",
    "length_biased"     = "Length-biased",
    "expression_biased" = "Expression-biased",
    "transcript_level"  = "Transcript-level",
    "id_corruption"     = "ID corruption"
  )

  plot_df$mode_label <- mode_labels[plot_df$mode]

  # Panel A: Global recovery recall
  p_recall <- ggplot(
    plot_df,
    aes(x = level_pct, y = mean_global_recall,
        color = method_label, fill = method_label, group = method_label)
  ) +
    geom_ribbon(
      aes(ymin = mean_global_recall - sd_global_recall,
          ymax = mean_global_recall + sd_global_recall),
      alpha = 0.15, color = NA
    ) +
    geom_line(linewidth = 1.0) +
    geom_point(size = 2.0) +
    facet_wrap(~ mode_label, ncol = 3) +
    scale_color_manual(values = method_colors) +
    scale_fill_manual(values = method_colors) +
    labs(
      title = "Global Recovery Recall by Degradation Mode",
      subtitle = "Fraction of full-annotation Tier A genes recovered",
      x = "Degradation level (%)",
      y = "Global recall",
      color = "Method",
      fill  = "Method"
    ) +
    theme_bw(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      strip.text      = element_text(size = 10, face = "bold"),
      panel.grid.minor = element_blank()
    )

  ggsave(output_fig, plot = p_recall, width = 12, height = 9, dpi = 150)
  cat(sprintf("  Wrote figure to %s\n", output_fig))

  # Panel B: Detectable recall (separate figure for clarity)
  p_detectable <- ggplot(
    plot_df,
    aes(x = level_pct, y = mean_detectable_recall,
        color = method_label, fill = method_label, group = method_label)
  ) +
    geom_ribbon(
      aes(ymin = mean_detectable_recall - sd_detectable_recall,
          ymax = mean_detectable_recall + sd_detectable_recall),
      alpha = 0.15, color = NA
    ) +
    geom_line(linewidth = 1.0) +
    geom_point(size = 2.0) +
    facet_wrap(~ mode_label, ncol = 3) +
    scale_color_manual(values = method_colors) +
    scale_fill_manual(values = method_colors) +
    labs(
      title = "Detectable Recall by Degradation Mode",
      subtitle = "Fraction of reference genes recovered (among genes still in degraded universe)",
      x = "Degradation level (%)",
      y = "Detectable recall",
      color = "Method",
      fill  = "Method"
    ) +
    theme_bw(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      strip.text      = element_text(size = 10, face = "bold"),
      panel.grid.minor = element_blank()
    )

  det_fig <- file.path(dirname(output_fig), "degradation_curve_detectable.png")
  ggsave(det_fig, plot = p_detectable, width = 12, height = 9, dpi = 150)
  cat(sprintf("  Wrote figure to %s\n", det_fig))

  # Panel C: Precision
  p_precision <- ggplot(
    plot_df,
    aes(x = level_pct, y = mean_precision,
        color = method_label, fill = method_label, group = method_label)
  ) +
    geom_ribbon(
      aes(ymin = mean_precision - sd_precision,
          ymax = mean_precision + sd_precision),
      alpha = 0.15, color = NA
    ) +
    geom_line(linewidth = 1.0) +
    geom_point(size = 2.0) +
    facet_wrap(~ mode_label, ncol = 3) +
    scale_color_manual(values = method_colors) +
    scale_fill_manual(values = method_colors) +
    labs(
      title = "Precision by Degradation Mode",
      subtitle = "Fraction of detected DEGs that are in the reference set",
      x = "Degradation level (%)",
      y = "Reference-consensus precision",
      color = "Method",
      fill  = "Method"
    ) +
    theme_bw(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      strip.text      = element_text(size = 10, face = "bold"),
      panel.grid.minor = element_blank()
    )

  prec_fig <- file.path(dirname(output_fig), "degradation_curve_precision.png")
  ggsave(prec_fig, plot = p_precision, width = 12, height = 9, dpi = 150)
  cat(sprintf("  Wrote figure to %s\n", prec_fig))

  # Combined multi-panel figure
  cat("\n── Generating combined multi-panel figure ──\n")

  # Use only consensus + featurecounts + salmon for cleaner visualization
  key_methods <- c("consensus", "featurecounts", "salmon")
  plot_sub <- plot_df[plot_df$method %in% key_methods, ]

  p1 <- ggplot(
    plot_sub,
    aes(x = level_pct, y = mean_global_recall,
        color = method_label, group = method_label)
  ) +
    geom_line(linewidth = 0.9) + geom_point(size = 1.8) +
    facet_wrap(~ mode_label, ncol = 5) +
    scale_color_manual(values = method_colors) +
    labs(title = "A: Global Recovery Recall", x = "Degradation (%)", y = "Recall") +
    theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "bottom", strip.text = element_text(size = 9, face = "bold"))

  p2 <- ggplot(
    plot_sub,
    aes(x = level_pct, y = mean_detectable_recall,
        color = method_label, group = method_label)
  ) +
    geom_line(linewidth = 0.9) + geom_point(size = 1.8) +
    facet_wrap(~ mode_label, ncol = 5) +
    scale_color_manual(values = method_colors) +
    labs(title = "B: Detectable Recall", x = "Degradation (%)", y = "Recall") +
    theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "bottom", strip.text = element_text(size = 9, face = "bold"))

  p3 <- ggplot(
    plot_sub,
    aes(x = level_pct, y = mean_precision,
        color = method_label, group = method_label)
  ) +
    geom_line(linewidth = 0.9) + geom_point(size = 1.8) +
    facet_wrap(~ mode_label, ncol = 5) +
    scale_color_manual(values = method_colors) +
    labs(title = "C: Precision", x = "Degradation (%)", y = "Reference-consensus precision") +
    theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "bottom", strip.text = element_text(size = 9, face = "bold"))

  # Combine using patchwork if available, else save separately and note
  combined_fig <- file.path(dirname(output_fig), "degradation_curve_combined.png")
  if (requireNamespace("patchwork", quietly = TRUE)) {
    library(patchwork)
    combo <- (p1 / p2 / p3) + patchwork::plot_layout(guides = "collect") &
      ggplot2::theme(legend.position = "bottom")
    ggsave(combined_fig, plot = combo, width = 14, height = 12, dpi = 150)
    cat(sprintf("  Wrote combined figure to %s\n", combined_fig))
  } else {
    ggsave(combined_fig, plot = p1, width = 14, height = 5, dpi = 150)
    cat(sprintf("  Wrote combined figure (panel A only) to %s\n", combined_fig))
    cat("  Note: install 'patchwork' package for full multi-panel combination.\n")
  }

  cat("\n========== Benchmark complete ==========\n")
}

main()
# nolint end
