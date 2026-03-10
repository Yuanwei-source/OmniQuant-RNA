# workflow/scripts/run_consensus_dea.R
# Consensus DEA across quantifiers using RRA as primary evidence and CCT as secondary evidence.

log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(RobustRankAggreg)
  library(UpSetR)
})

invisible(utils::globalVariables(c(
  "gene_id_standard", "gene_name", "logFC", "P.Value", "adj.P.Val",
  "included_in_main", "method", "contrast"
)))

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

as_flag <- function(x, default = FALSE) {
  if (is.null(x) || length(x) == 0 || is.na(x[1])) {
    return(default)
  }
  if (is.logical(x)) {
    return(x)
  }
  tolower(as.character(x)) %in% c("true", "1", "yes")
}

clip_probabilities <- function(p, eps) {
  pmin(pmax(p, eps), 1 - eps)
}

cct_pvalue <- function(p_values, eps) {
  if (length(p_values) == 0) {
    return(1)
  }
  p_values <- clip_probabilities(as.numeric(p_values), eps)
  t_stat <- mean(tan((0.5 - p_values) * pi))
  combined <- 0.5 - atan(t_stat) / pi
  pmin(pmax(combined, eps), 1)
}

build_rank_list <- function(df, direction) {
  ranked <- df %>%
    filter(!is.na(.data$gene_id_standard), !is.na(.data$P.Value), !is.na(.data$logFC))

  if (direction == "up") {
    ranked <- ranked %>% filter(.data$logFC > 0) %>% arrange(.data$P.Value, desc(.data$logFC), .data$gene_id_standard)
  } else {
    ranked <- ranked %>% filter(.data$logFC < 0) %>% arrange(.data$P.Value, .data$logFC, .data$gene_id_standard)
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
    exact = length(rank_lists) < 10
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

compute_cct_scores <- function(logfc_matrix, p_matrix, quantifiers, universe, direction, eps) {
  cct_p <- vapply(
    seq_len(nrow(logfc_matrix)),
    function(i) {
      logfc_values <- logfc_matrix[i, ]
      p_values <- p_matrix[i, ]
      keep <- !is.na(logfc_values) & !is.na(p_values)
      if (direction == "up") {
        keep <- keep & logfc_values > 0
      } else {
        keep <- keep & logfc_values < 0
      }
      cct_pvalue(p_values[keep], eps = eps)
    },
    numeric(1)
  )

  tibble(
    gene_id_standard = universe,
    p = cct_p,
    fdr = p.adjust(cct_p, method = "BH")
  )
}

coalesce_columns <- function(df, cols) {
  if (length(cols) == 0) {
    return(rep(NA_character_, nrow(df)))
  }

  current <- df[[cols[1]]]
  if (length(cols) == 1) {
    return(current)
  }

  for (col in cols[-1]) {
    current <- dplyr::coalesce(current, df[[col]])
  }
  current
}

collapse_names <- function(values, labels) {
  labels[!is.na(values)]
}

choose_direction <- function(up_support, down_support, rra_up_p, rra_down_p) {
  if (up_support == 0 && down_support == 0) {
    return("none")
  }
  if (up_support > down_support) {
    return("up")
  }
  if (down_support > up_support) {
    return("down")
  }
  if (rra_up_p < rra_down_p) {
    return("up")
  }
  if (rra_down_p < rra_up_p) {
    return("down")
  }
  "mixed"
}

compute_consensus_logfc <- function(values, direction) {
  values <- values[!is.na(values)]
  if (length(values) == 0) {
    return(NA_real_)
  }
  if (direction == "up") {
    values <- values[values > 0]
  } else if (direction == "down") {
    values <- values[values < 0]
  }
  if (length(values) == 0) {
    return(NA_real_)
  }
  median(values)
}

compute_logfc_cv <- function(values, direction) {
  values <- values[!is.na(values)]
  if (length(values) < 2) {
    return(NA_real_)
  }

  if (direction == "up") {
    values <- values[values > 0]
  } else if (direction == "down") {
    values <- abs(values[values < 0])
  } else {
    values <- abs(values)
  }

  if (length(values) < 2) {
    return(NA_real_)
  }

  mean_value <- mean(values)
  if (is.na(mean_value) || mean_value == 0) {
    return(NA_real_)
  }

  sd(values) / mean_value
}

assign_tier <- function(support_n, sign_consistency_n, best_rra_fdr, best_cct_fdr, logfc_cv, tiers) {
  if (
    support_n >= tiers$tier_a$min_support &&
    sign_consistency_n >= tiers$tier_a$min_sign_consistency &&
    best_rra_fdr <= tiers$tier_a$max_rra_fdr &&
    best_cct_fdr <= tiers$tier_a$max_cct_fdr &&
    (is.na(logfc_cv) || logfc_cv <= tiers$tier_a$max_logfc_cv)
  ) {
    return("Tier_A")
  }

  if (
    support_n >= tiers$tier_b$min_support &&
    sign_consistency_n >= tiers$tier_b$min_sign_consistency &&
    best_rra_fdr <= tiers$tier_b$max_rra_fdr &&
    best_cct_fdr <= tiers$tier_b$max_cct_fdr &&
    (is.na(logfc_cv) || logfc_cv <= tiers$tier_b$max_logfc_cv)
  ) {
    return("Tier_B")
  }

  if (
    support_n >= tiers$tier_c$min_support &&
    sign_consistency_n >= tiers$tier_c$min_sign_consistency &&
    best_rra_fdr <= tiers$tier_c$max_rra_fdr &&
    best_cct_fdr <= tiers$tier_c$max_cct_fdr &&
    (is.na(logfc_cv) || logfc_cv <= tiers$tier_c$max_logfc_cv)
  ) {
    return("Tier_C")
  }

  "unclassified"
}

list_failures <- function(support_n, sign_consistency_n, best_rra_fdr, best_cct_fdr, logfc_cv, tier_cfg) {
  failures <- c()
  if (support_n < tier_cfg$min_support) failures <- c(failures, "support")
  if (sign_consistency_n < tier_cfg$min_sign_consistency) failures <- c(failures, "sign_consistency")
  if (best_rra_fdr > tier_cfg$max_rra_fdr) failures <- c(failures, "rra_fdr")
  if (best_cct_fdr > tier_cfg$max_cct_fdr) failures <- c(failures, "cct_fdr")
  if (!is.na(logfc_cv) && logfc_cv > tier_cfg$max_logfc_cv) failures <- c(failures, "logFC_CV")
  if (length(failures) == 0) "pass" else paste(failures, collapse = ";")
}

classify_tier_with_override <- function(support_n, sign_consistency_n, best_rra_fdr, best_cct_fdr, logfc_cv, tiers,
                                        tier_b_support = NULL, ignore_cv = FALSE, ignore_cct = FALSE, ignore_rra = FALSE) {
  tiers_local <- tiers
  if (!is.null(tier_b_support)) {
    tiers_local$tier_b$min_support <- tier_b_support
  }
  if (ignore_cv) {
    tiers_local$tier_a$max_logfc_cv <- Inf
    tiers_local$tier_b$max_logfc_cv <- Inf
    tiers_local$tier_c$max_logfc_cv <- Inf
  }
  if (ignore_cct) {
    tiers_local$tier_a$max_cct_fdr <- Inf
    tiers_local$tier_b$max_cct_fdr <- Inf
    tiers_local$tier_c$max_cct_fdr <- Inf
  }
  if (ignore_rra) {
    tiers_local$tier_a$max_rra_fdr <- Inf
    tiers_local$tier_b$max_rra_fdr <- Inf
    tiers_local$tier_c$max_rra_fdr <- Inf
  }
  assign_tier(support_n, sign_consistency_n, best_rra_fdr, best_cct_fdr, logfc_cv, tiers_local)
}

load_quantifier_result <- function(path, quantifier, method, contrast) {
  read_csv(path, show_col_types = FALSE) %>%
    filter(
      .data$method == method,
      .data$contrast == contrast,
      as_flag(.data$included_in_main, default = TRUE)
    ) %>%
    mutate(
      gene_id_standard = coalesce(.data$gene_id_standard, .data$gene_id),
      gene_name = coalesce(.data$gene_name, .data$gene_id_standard),
      quantifier = quantifier
    ) %>%
    filter(!is.na(.data$gene_id_standard), .data$gene_id_standard != "") %>%
    arrange(.data$adj.P.Val, .data$P.Value, desc(abs(.data$logFC))) %>%
    distinct(.data$gene_id_standard, .keep_all = TRUE) %>%
    select("gene_id_standard", "gene_name", "logFC", "P.Value", "adj.P.Val", "quantifier")
}

resolve_quantifier_result_path <- function(directory_path, method, contrast) {
  file_path <- file.path(directory_path, paste0(method, ".", contrast, ".csv"))
  if (!file.exists(file_path)) {
    stop("Missing DEA contrast file: ", file_path)
  }
  file_path
}

plot_scatter <- function(consensus_df, output_path, contrast) {
  x_col <- "logFC__featurecounts"
  y_col <- "logFC__salmon"

  if (!(x_col %in% colnames(consensus_df) && y_col %in% colnames(consensus_df))) {
    pdf(output_path, width = 7, height = 6)
    plot.new()
    text(0.5, 0.5, "FeatureCounts or Salmon logFC columns are missing")
    dev.off()
    return(invisible(NULL))
  }

  plot_df <- consensus_df %>%
    transmute(
      gene_id_standard = .data$gene_id_standard,
      featurecounts = .data[[x_col]],
      salmon = .data[[y_col]],
      tier = .data$tier
    ) %>%
    filter(!is.na(.data$featurecounts), !is.na(.data$salmon))

  if (nrow(plot_df) == 0) {
    pdf(output_path, width = 7, height = 6)
    plot.new()
    text(0.5, 0.5, "No overlapping logFC values for Salmon vs featureCounts")
    dev.off()
    return(invisible(NULL))
  }

  correlation <- suppressWarnings(cor(plot_df$featurecounts, plot_df$salmon, use = "complete.obs"))
  p <- ggplot(plot_df, aes(x = .data$featurecounts, y = .data$salmon, color = .data$tier)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(alpha = 0.75, size = 1.8) +
    labs(
      title = paste0("logFC consistency: Salmon vs featureCounts (", contrast, ")"),
      subtitle = paste0("Pearson r = ", round(correlation, 3), "; n = ", nrow(plot_df)),
      x = "featureCounts logFC",
      y = "Salmon logFC",
      color = "Tier"
    ) +
    theme_bw(base_size = 11)

  ggsave(output_path, plot = p, width = 7, height = 6)
}

plot_consensus_volcano <- function(consensus_df, output_path, contrast, p_clip) {
  plot_df <- consensus_df %>%
    mutate(
      volcano_y = -log10(pmax(.data$best_rra_fdr, p_clip)),
      tier = factor(.data$tier, levels = c("Tier_A", "Tier_B", "Tier_C", "unclassified"))
    )

  p <- ggplot(plot_df, aes(x = .data$consensus_logFC, y = .data$volcano_y, color = .data$tier)) +
    geom_point(alpha = 0.75, size = 1.6) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    labs(
      title = paste0("Consensus volcano: ", contrast),
      x = "consensus_logFC",
      y = "-log10(best_rra_fdr)",
      color = "Tier"
    ) +
    theme_bw(base_size = 11)

  ggsave(output_path, plot = p, width = 7, height = 6)
}

plot_significance_upset <- function(membership_df, quantifiers, output_path, contrast) {
  sets <- lapply(quantifiers, function(quantifier) membership_df$gene_id_standard[membership_df[[quantifier]]])
  names(sets) <- quantifiers
  non_empty_sets <- sets[lengths(sets) > 0]

  pdf(output_path, width = 9, height = 6)
  on.exit(dev.off(), add = TRUE)

  if (length(non_empty_sets) == 0) {
    plot.new()
    text(0.5, 0.5, paste0("No significant genes for UpSet: ", contrast))
    return(invisible(NULL))
  }

  upset(
    fromList(non_empty_sets),
    nsets = length(non_empty_sets),
    nintersects = NA,
    order.by = c("freq", "degree"),
    mainbar.y.label = "Significant gene intersections",
    sets.x.label = "Significant genes per quantifier",
    text.scale = 1.1
  )
}

build_sensitivity_table <- function(consensus_df, tiers) {
  tibble(
    scenario = c(
      "current_tier_b",
      "tier_b_support2",
      "tier_b_support4",
      "tier_b_without_cv",
      "tier_b_without_cct",
      "tier_b_without_rra"
    ),
    genes = c(
      sum(consensus_df$tier == "Tier_B"),
      sum(vapply(seq_len(nrow(consensus_df)), function(i) classify_tier_with_override(consensus_df$support_n[i], consensus_df$sign_consistency_n[i], consensus_df$best_rra_fdr[i], consensus_df$best_cct_fdr[i], consensus_df$logFC_CV[i], tiers, tier_b_support = 2) == "Tier_B", logical(1))),
      sum(vapply(seq_len(nrow(consensus_df)), function(i) classify_tier_with_override(consensus_df$support_n[i], consensus_df$sign_consistency_n[i], consensus_df$best_rra_fdr[i], consensus_df$best_cct_fdr[i], consensus_df$logFC_CV[i], tiers, tier_b_support = 4) == "Tier_B", logical(1))),
      sum(vapply(seq_len(nrow(consensus_df)), function(i) classify_tier_with_override(consensus_df$support_n[i], consensus_df$sign_consistency_n[i], consensus_df$best_rra_fdr[i], consensus_df$best_cct_fdr[i], consensus_df$logFC_CV[i], tiers, ignore_cv = TRUE) == "Tier_B", logical(1))),
      sum(vapply(seq_len(nrow(consensus_df)), function(i) classify_tier_with_override(consensus_df$support_n[i], consensus_df$sign_consistency_n[i], consensus_df$best_rra_fdr[i], consensus_df$best_cct_fdr[i], consensus_df$logFC_CV[i], tiers, ignore_cct = TRUE) == "Tier_B", logical(1))),
      sum(vapply(seq_len(nrow(consensus_df)), function(i) classify_tier_with_override(consensus_df$support_n[i], consensus_df$sign_consistency_n[i], consensus_df$best_rra_fdr[i], consensus_df$best_cct_fdr[i], consensus_df$logFC_CV[i], tiers, ignore_rra = TRUE) == "Tier_B", logical(1)))
    )
  )
}

consensus_config <- snakemake@config[["consensus"]]
config_dea <- snakemake@config[["dea"]]
contrast <- snakemake@wildcards[["contrast"]]
method <- consensus_config$methods[[1]] %||% "deseq2"
configured_quantifiers <- consensus_config$quantifiers %||% c("featurecounts", "stringtie", "salmon", "kallisto")
p_clip <- as.numeric(consensus_config$p_clip %||% 1e-16)
tiers <- consensus_config$tiers
fdr_threshold <- as.numeric(config_dea$fdr_threshold %||% 0.05)
lfc_threshold <- as.numeric(config_dea$lfc_threshold %||% 1.0)

raw_input_paths <- unique(unname(unlist(snakemake@input)))
raw_input_paths <- raw_input_paths[!is.na(raw_input_paths)]
if (length(raw_input_paths) < length(configured_quantifiers)) {
  stop("Consensus inputs are incomplete after extracting Snakemake paths.")
}
raw_input_paths <- raw_input_paths[seq_along(configured_quantifiers)]
quantifiers <- configured_quantifiers
input_paths <- as.list(raw_input_paths)
names(input_paths) <- quantifiers
output_table <- snakemake@output[["table"]]
output_summary <- snakemake@output[["summary"]]
output_diagnostics <- snakemake@output[["diagnostics"]]
output_membership <- snakemake@output[["membership"]]
output_sensitivity <- snakemake@output[["sensitivity"]]
output_scatter <- snakemake@output[["scatter"]]
output_volcano <- snakemake@output[["volcano"]]
output_upset <- snakemake@output[["upset"]]

dir.create(dirname(output_table), recursive = TRUE, showWarnings = FALSE)

cat("Loading consensus inputs for contrast:", contrast, "\n")
cat("Consensus input quantifiers:", paste(names(input_paths), collapse = ", "), "\n")
data_list <- lapply(names(input_paths), function(quantifier) {
  contrast_path <- resolve_quantifier_result_path(input_paths[[quantifier]], method, contrast)
  load_quantifier_result(contrast_path, quantifier, method, contrast)
})
names(data_list) <- names(input_paths)

universe <- sort(unique(unlist(lapply(data_list, function(df) df$gene_id_standard))))
if (length(universe) == 0) {
  stop("No genes available for consensus after main-namespace filtering.")
}

cat("Consensus universe size N =", length(universe), "\n")

consensus_df <- tibble(gene_id_standard = universe)
for (quantifier in names(data_list)) {
  df <- data_list[[quantifier]] %>%
    transmute(
      gene_id_standard = .data$gene_id_standard,
      !!paste0("gene_name__", quantifier) := .data$gene_name,
      !!paste0("logFC__", quantifier) := .data$logFC,
      !!paste0("P.Value__", quantifier) := .data$P.Value,
      !!paste0("adj.P.Val__", quantifier) := .data$adj.P.Val
    )
  consensus_df <- consensus_df %>% left_join(df, by = "gene_id_standard")
}

gene_name_cols <- grep("^gene_name__", colnames(consensus_df), value = TRUE)
logfc_cols <- paste0("logFC__", names(data_list))
p_cols <- paste0("P.Value__", names(data_list))
adj_cols <- paste0("adj.P.Val__", names(data_list))

consensus_df$gene_name <- coalesce_columns(consensus_df, gene_name_cols)

logfc_matrix <- as.matrix(consensus_df[, logfc_cols, drop = FALSE])
p_matrix <- as.matrix(consensus_df[, p_cols, drop = FALSE])
adj_matrix <- as.matrix(consensus_df[, adj_cols, drop = FALSE])
storage.mode(logfc_matrix) <- "numeric"
storage.mode(p_matrix) <- "numeric"
storage.mode(adj_matrix) <- "numeric"

sig_matrix <- !is.na(logfc_matrix) & !is.na(adj_matrix) & (adj_matrix <= fdr_threshold) & (abs(logfc_matrix) >= lfc_threshold)

support_n <- rowSums(sig_matrix, na.rm = TRUE)
up_support_n <- rowSums(sig_matrix & logfc_matrix > 0, na.rm = TRUE)
down_support_n <- rowSums(sig_matrix & logfc_matrix < 0, na.rm = TRUE)
sign_consistency_n <- pmax(up_support_n, down_support_n)

rra_up <- compute_rra_scores(data_list, universe, direction = "up") %>% rename(rra_up_p = "p", rra_up_fdr = "fdr")
rra_down <- compute_rra_scores(data_list, universe, direction = "down") %>% rename(rra_down_p = "p", rra_down_fdr = "fdr")
cct_up <- compute_cct_scores(logfc_matrix, p_matrix, names(data_list), universe, direction = "up", eps = p_clip) %>% rename(cct_up_p = "p", cct_up_fdr = "fdr")
cct_down <- compute_cct_scores(logfc_matrix, p_matrix, names(data_list), universe, direction = "down", eps = p_clip) %>% rename(cct_down_p = "p", cct_down_fdr = "fdr")

consensus_df <- consensus_df %>%
  left_join(rra_up, by = "gene_id_standard") %>%
  left_join(rra_down, by = "gene_id_standard") %>%
  left_join(cct_up, by = "gene_id_standard") %>%
  left_join(cct_down, by = "gene_id_standard")

consensus_direction <- vapply(
  seq_len(nrow(consensus_df)),
  function(i) choose_direction(up_support_n[i], down_support_n[i], consensus_df$rra_up_p[i], consensus_df$rra_down_p[i]),
  character(1)
)

supporting_quantifiers <- vapply(
  seq_len(nrow(sig_matrix)),
  function(i) paste(names(data_list)[sig_matrix[i, ]], collapse = ";"),
  character(1)
)

consistent_quantifiers <- vapply(
  seq_len(nrow(sig_matrix)),
  function(i) {
    direction <- consensus_direction[i]
    if (direction == "up") {
      keep <- sig_matrix[i, ] & logfc_matrix[i, ] > 0
    } else if (direction == "down") {
      keep <- sig_matrix[i, ] & logfc_matrix[i, ] < 0
    } else {
      keep <- sig_matrix[i, ]
    }
    paste(names(data_list)[keep], collapse = ";")
  },
  character(1)
)

consensus_logfc <- vapply(
  seq_len(nrow(logfc_matrix)),
  function(i) compute_consensus_logfc(logfc_matrix[i, ], consensus_direction[i]),
  numeric(1)
)

logfc_mean <- rowMeans(logfc_matrix, na.rm = TRUE)
logfc_mean[is.nan(logfc_mean)] <- NA_real_

logfc_median <- vapply(
  seq_len(nrow(logfc_matrix)),
  function(i) {
    values <- logfc_matrix[i, ]
    values <- values[!is.na(values)]
    if (length(values) == 0) NA_real_ else median(values)
  },
  numeric(1)
)

logfc_sd <- vapply(
  seq_len(nrow(logfc_matrix)),
  function(i) {
    values <- logfc_matrix[i, ]
    values <- values[!is.na(values)]
    if (length(values) < 2) NA_real_ else sd(values)
  },
  numeric(1)
)

logfc_cv <- vapply(
  seq_len(nrow(logfc_matrix)),
  function(i) compute_logfc_cv(logfc_matrix[i, ], consensus_direction[i]),
  numeric(1)
)

best_rra_p <- ifelse(consensus_direction == "down", consensus_df$rra_down_p, consensus_df$rra_up_p)
best_rra_fdr <- ifelse(consensus_direction == "down", consensus_df$rra_down_fdr, consensus_df$rra_up_fdr)
best_cct_p <- ifelse(consensus_direction == "down", consensus_df$cct_down_p, consensus_df$cct_up_p)
best_cct_fdr <- ifelse(consensus_direction == "down", consensus_df$cct_down_fdr, consensus_df$cct_up_fdr)

best_rra_p[consensus_direction == "mixed"] <- pmin(consensus_df$rra_up_p[consensus_direction == "mixed"], consensus_df$rra_down_p[consensus_direction == "mixed"])
best_rra_fdr[consensus_direction == "mixed"] <- pmin(consensus_df$rra_up_fdr[consensus_direction == "mixed"], consensus_df$rra_down_fdr[consensus_direction == "mixed"])
best_cct_p[consensus_direction == "mixed"] <- pmin(consensus_df$cct_up_p[consensus_direction == "mixed"], consensus_df$cct_down_p[consensus_direction == "mixed"])
best_cct_fdr[consensus_direction == "mixed"] <- pmin(consensus_df$cct_up_fdr[consensus_direction == "mixed"], consensus_df$cct_down_fdr[consensus_direction == "mixed"])

tier <- vapply(
  seq_len(nrow(consensus_df)),
  function(i) {
    assign_tier(
      support_n = support_n[i],
      sign_consistency_n = sign_consistency_n[i],
      best_rra_fdr = best_rra_fdr[i],
      best_cct_fdr = best_cct_fdr[i],
      logfc_cv = logfc_cv[i],
      tiers = tiers
    )
  },
  character(1)
)

tier_blocker_a <- vapply(
  seq_len(nrow(consensus_df)),
  function(i) list_failures(support_n[i], sign_consistency_n[i], best_rra_fdr[i], best_cct_fdr[i], logfc_cv[i], tiers$tier_a),
  character(1)
)

tier_blocker_b <- vapply(
  seq_len(nrow(consensus_df)),
  function(i) list_failures(support_n[i], sign_consistency_n[i], best_rra_fdr[i], best_cct_fdr[i], logfc_cv[i], tiers$tier_b),
  character(1)
)

tier_blocker_c <- vapply(
  seq_len(nrow(consensus_df)),
  function(i) list_failures(support_n[i], sign_consistency_n[i], best_rra_fdr[i], best_cct_fdr[i], logfc_cv[i], tiers$tier_c),
  character(1)
)

membership_df <- tibble(gene_id_standard = universe, gene_name = consensus_df$gene_name)
for (idx in seq_along(quantifiers)) {
  membership_df[[quantifiers[idx]]] <- sig_matrix[, idx]
}
membership_df <- membership_df %>%
  mutate(
    support_n = support_n,
    up_support_n = up_support_n,
    down_support_n = down_support_n,
    membership_pattern = vapply(seq_len(nrow(sig_matrix)), function(i) paste(quantifiers[sig_matrix[i, ]], collapse = ";"), character(1))
  )

consensus_df <- consensus_df %>%
  mutate(
    contrast = contrast,
    method = method,
    universe_n = length(universe),
    support_n = support_n,
    up_support_n = up_support_n,
    down_support_n = down_support_n,
    sign_consistency_n = sign_consistency_n,
    consensus_direction = consensus_direction,
    supporting_quantifiers = supporting_quantifiers,
    consistent_quantifiers = consistent_quantifiers,
    consensus_logFC = consensus_logfc,
    logFC_mean = logfc_mean,
    logFC_median = logfc_median,
    logFC_sd = logfc_sd,
    logFC_CV = logfc_cv,
    best_rra_p = best_rra_p,
    best_rra_fdr = best_rra_fdr,
    best_cct_p = best_cct_p,
    best_cct_fdr = best_cct_fdr,
    tier = tier,
    tier_blocker_a = tier_blocker_a,
    tier_blocker_b = tier_blocker_b,
    tier_blocker_c = tier_blocker_c
  ) %>%
  select(
    "gene_id_standard",
    "gene_name",
    "contrast",
    "method",
    "universe_n",
    "support_n",
    "up_support_n",
    "down_support_n",
    "sign_consistency_n",
    "consensus_direction",
    "supporting_quantifiers",
    "consistent_quantifiers",
    "consensus_logFC",
    "logFC_mean",
    "logFC_median",
    "logFC_sd",
    "logFC_CV",
    "rra_up_p",
    "rra_up_fdr",
    "rra_down_p",
    "rra_down_fdr",
    "cct_up_p",
    "cct_up_fdr",
    "cct_down_p",
    "cct_down_fdr",
    "best_rra_p",
    "best_rra_fdr",
    "best_cct_p",
    "best_cct_fdr",
    "tier",
    "tier_blocker_a",
    "tier_blocker_b",
    "tier_blocker_c",
    all_of(logfc_cols),
    all_of(p_cols),
    all_of(adj_cols)
  ) %>%
  arrange(.data$best_rra_fdr, .data$best_cct_fdr, desc(.data$sign_consistency_n), desc(abs(.data$consensus_logFC)))

diagnostics_df <- consensus_df %>%
  transmute(
    gene_id_standard = .data$gene_id_standard,
    gene_name = .data$gene_name,
    tier = .data$tier,
    support_n = .data$support_n,
    sign_consistency_n = .data$sign_consistency_n,
    consensus_direction = .data$consensus_direction,
    best_rra_fdr = .data$best_rra_fdr,
    best_cct_fdr = .data$best_cct_fdr,
    logFC_CV = .data$logFC_CV,
    tier_blocker_a = .data$tier_blocker_a,
    tier_blocker_b = .data$tier_blocker_b,
    tier_blocker_c = .data$tier_blocker_c
  )

sensitivity_df <- build_sensitivity_table(consensus_df, tiers)

summary_df <- tibble(
  metric = c(
    "contrast",
    "method",
    "quantifiers",
    "universe_n",
    "significance_fdr_threshold",
    "significance_lfc_threshold",
    "tier_a_n",
    "tier_b_n",
    "tier_c_n",
    "unclassified_n",
    "support2_tier_b_n",
    "tier_b_without_cv_n",
    "tier_b_without_cct_n",
    "tier_b_without_rra_n"
  ),
  value = c(
    contrast,
    method,
    paste(names(data_list), collapse = ";"),
    length(universe),
    fdr_threshold,
    lfc_threshold,
    sum(consensus_df$tier == "Tier_A"),
    sum(consensus_df$tier == "Tier_B"),
    sum(consensus_df$tier == "Tier_C"),
    sum(consensus_df$tier == "unclassified"),
    sensitivity_df$genes[sensitivity_df$scenario == "tier_b_support2"],
    sensitivity_df$genes[sensitivity_df$scenario == "tier_b_without_cv"],
    sensitivity_df$genes[sensitivity_df$scenario == "tier_b_without_cct"],
    sensitivity_df$genes[sensitivity_df$scenario == "tier_b_without_rra"]
  )
)

write_tsv(consensus_df, output_table)
write_tsv(summary_df, output_summary)
write_tsv(diagnostics_df, output_diagnostics)
write_tsv(membership_df, output_membership)
write_tsv(sensitivity_df, output_sensitivity)
plot_scatter(consensus_df, output_scatter, contrast)
plot_consensus_volcano(consensus_df, output_volcano, contrast, p_clip)
plot_significance_upset(membership_df, quantifiers, output_upset, contrast)

cat("Consensus DEA completed successfully.\n")
sink()