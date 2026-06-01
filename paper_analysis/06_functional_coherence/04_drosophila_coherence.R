#!/usr/bin/env Rscript
###########################################################################
# Drosophila GO Enrichment: Tier A vs near-Tier-C vs random control
# Uses org.Dm.eg.db + clusterProfiler for GO BP enrichment
# Size-matched random control (1,000 iterations) via fast hypergeometric
###########################################################################
suppressPackageStartupMessages({
  library(org.Dm.eg.db)
  library(clusterProfiler)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(GO.db)
})

source("../theme_nature.R")

# ── Paths ──────────────────────────────────────────────────────────────
BASE_DIR    <- Sys.getenv("PROJECT_ROOT", normalizePath("."))
EXP_DIR     <- file.path(BASE_DIR, "experiments/drosophila_enrichment")
ON_TSV      <- file.path(BASE_DIR, "benchmark_results/drosophila_wolbachia/runs/2026-05-27_drosophila_wolbachia_decontam-on/07.consensus_expression/Wolbachia_infected_vs_Wolbachia_free/consensus_results.tsv")
OFF_TSV     <- file.path(BASE_DIR, "benchmark_results/drosophila_wolbachia/runs/2026-05-27_drosophila_wolbachia_decontam-off/07.consensus_expression/Wolbachia_infected_vs_Wolbachia_free/consensus_results.tsv")

# Output dirs
GENESETS_DIR <- file.path(EXP_DIR, "genesets")
TABLES_DIR   <- file.path(EXP_DIR, "results/tables")
FIGURES_DIR  <- file.path(EXP_DIR, "results/figures")

# Enrichment params
ONT          <- "BP"
PVAL_CUTOFF  <- 0.05
QVALUE_CUTOFF <- 0.2
MIN_GS_SIZE  <- 5
MAX_GS_SIZE  <- 500
N_ITER       <- 1000  # random control iterations
FDR_THRESH   <- 0.05  # significance threshold for counting terms

cat("═══ Drosophila GO Enrichment Analysis ═══\n\n")

# ═══════════════════════════════════════════════════════════════════════
# STEP 1: Load and prepare data
# ═══════════════════════════════════════════════════════════════════════
cat("── Step 1: Loading consensus data ──\n")
on  <- read_tsv(ON_TSV,  show_col_types = FALSE)
off <- read_tsv(OFF_TSV, show_col_types = FALSE)

cat(sprintf("  decontam-ON: %d genes\n", nrow(on)))
cat(sprintf("  decontam-OFF: %d genes\n", nrow(off)))

# ═══════════════════════════════════════════════════════════════════════
# STEP 2: Define gene sets
# ═══════════════════════════════════════════════════════════════════════
cat("\n── Step 2: Defining gene sets ──\n")

# Helper: strip "gene:" prefix → FBgn ID
strip_prefix <- function(x) sub("^gene:", "", x)

# Tier A (decontam-ON)
tier_a_ids <- on %>%
  filter(tier == "Tier_A") %>%
  pull(gene_id_standard) %>%
  unique() %>%
  strip_prefix()
cat(sprintf("  Tier A: %d genes\n", length(tier_a_ids)))

# near-Tier-C: support_n >= 2 AND unclassified (NOT Tier A/B/C)
near_tc_ids <- on %>%
  filter(support_n >= 2, tier == "unclassified") %>%
  pull(gene_id_standard) %>%
  unique() %>%
  strip_prefix()
cat(sprintf("  near-Tier-C: %d genes\n", length(near_tc_ids)))

# OFF-only: Tier_A in OFF but NOT Tier_A in ON
off_tier_a <- off %>% filter(tier == "Tier_A") %>% 
  pull(gene_id_standard) %>% unique() %>% strip_prefix()
on_tier_a  <- on  %>% filter(tier == "Tier_A") %>% 
  pull(gene_id_standard) %>% unique() %>% strip_prefix()
off_only_ids <- setdiff(off_tier_a, on_tier_a)
cat(sprintf("  OFF-only: %d genes\n", length(off_only_ids)))

# ═══════════════════════════════════════════════════════════════════════
# STEP 3: Pre-compute GO annotation
# ═══════════════════════════════════════════════════════════════════════
cat("\n── Step 3: Pre-computing GO annotations ──\n")

all_maps <- AnnotationDbi::select(org.Dm.eg.db,
  keys    = keys(org.Dm.eg.db, "FLYBASE"),
  columns = c("GO", "ONTOLOGY"),
  keytype = "FLYBASE"
)
# Filter to BP only
go_maps <- all_maps[all_maps$ONTOLOGY == ONT & !is.na(all_maps$GO), c("GO", "FLYBASE")]
rm(all_maps); invisible(gc())

# Filter terms by size
term_counts <- table(go_maps$GO)
valid_terms <- names(term_counts)[term_counts >= MIN_GS_SIZE & term_counts <= MAX_GS_SIZE]
go_maps <- go_maps[go_maps$GO %in% valid_terms, ]
cat(sprintf("  GO BP terms (5-500 genes): %d\n", length(valid_terms)))
cat(sprintf("  Annotated genes: %d\n", length(unique(go_maps$FLYBASE))))

# Build TERM2GENE and TERM2NAME for enricher()
term2gene <- go_maps[, c("GO", "FLYBASE")]
term2name <- data.frame(
  GO   = valid_terms,
  Term = AnnotationDbi::mapIds(GO.db, keys = valid_terms, column = "TERM", keytype = "GOID"),
  stringsAsFactors = FALSE
)
# Fill any missing names
term2name$Term[is.na(term2name$Term)] <- term2name$GO[is.na(term2name$Term)]

# Build term-to-genes list for fast random control
term2genes_list <- split(go_maps$FLYBASE, go_maps$GO)
M_vec <- lengths(term2genes_list)            # number of genes per term
all_bg_genes <- unique(go_maps$FLYBASE)     # background gene universe
N_bg <- length(all_bg_genes)
cat(sprintf("  Background: %d genes, %d terms\n", N_bg, length(term2genes_list)))

# ═══════════════════════════════════════════════════════════════════════
# STEP 4: Filter gene sets to annotated genes + save gene lists
# ═══════════════════════════════════════════════════════════════════════
cat("\n── Step 4: Filtering gene sets ──\n")

tier_a_annot    <- intersect(tier_a_ids, all_bg_genes)
near_tc_annot   <- intersect(near_tc_ids, all_bg_genes)
off_only_annot  <- intersect(off_only_ids, all_bg_genes)

cat(sprintf("  Tier A annotated: %d / %d (%.1f%%)\n",
  length(tier_a_annot), length(tier_a_ids),
  100*length(tier_a_annot)/length(tier_a_ids)))
cat(sprintf("  near-Tier-C annotated: %d / %d (%.1f%%)\n",
  length(near_tc_annot), length(near_tc_ids),
  100*length(near_tc_annot)/length(near_tc_ids)))
cat(sprintf("  OFF-only annotated: %d / %d (%.1f%%)\n",
  length(off_only_annot), length(off_only_ids),
  100*length(off_only_annot)/length(off_only_ids)))

# Save gene lists
writeLines(tier_a_ids,   file.path(GENESETS_DIR, "tier_a.txt"))
writeLines(near_tc_ids,  file.path(GENESETS_DIR, "near_tier_c.txt"))
writeLines(off_only_ids, file.path(GENESETS_DIR, "off_only.txt"))
writeLines(all_bg_genes, file.path(GENESETS_DIR, "background.txt"))
cat("  Gene lists saved to genesets/\n")

# ═══════════════════════════════════════════════════════════════════════
# STEP 5: GO Enrichment for Tier A, near-Tier-C, OFF-only
# ═══════════════════════════════════════════════════════════════════════
cat("\n── Step 5: Running GO enrichment ──\n")

run_enrichment <- function(genes, label) {
  cat(sprintf("  %s (%d genes)...\n", label, length(genes)))
  if (length(genes) < 5) {
    warning(sprintf("  %s has < 5 annotated genes, skipping", label))
    return(NULL)
  }
  er <- enricher(
    gene          = genes,
    universe      = all_bg_genes,
    TERM2GENE     = term2gene,
    TERM2NAME     = term2name,
    pAdjustMethod = "BH",
    pvalueCutoff  = PVAL_CUTOFF,
    qvalueCutoff  = QVALUE_CUTOFF,
    minGSSize     = MIN_GS_SIZE,
    maxGSSize     = MAX_GS_SIZE
  )
  if (is.null(er)) {
    cat(sprintf("    No significant terms.\n"))
    return(data.frame())
  }
  res <- as.data.frame(er)
  cat(sprintf("    %d significant terms\n", nrow(res)))
  res
}

tier_a_res   <- run_enrichment(tier_a_annot,   "Tier A")
near_tc_res  <- run_enrichment(near_tc_annot,  "near-Tier-C")
off_only_res <- run_enrichment(off_only_annot, "OFF-only")

# ═══════════════════════════════════════════════════════════════════════
# STEP 6: Size-matched random control (1,000 iterations)
# Uses fast vectorized hypergeometric test
# ═══════════════════════════════════════════════════════════════════════
cat(sprintf("\n── Step 6: Random control (%d iterations) ──\n", N_ITER))

# Match Tier A size
n_target <- length(tier_a_annot)
cat(sprintf("  Sampling %d genes per iteration (matching Tier A)\n", n_target))

set.seed(20260529)

# Pre-allocate result vector
iter_sig_counts <- integer(N_ITER)

# Progress tracking
pb_step <- max(1, floor(N_ITER / 20))

t_start <- Sys.time()
for (i in seq_len(N_ITER)) {
  # Random sample from background
  sample_genes <- sample(all_bg_genes, n_target, replace = FALSE)
  
  # Compute overlap (k) for all terms
  k_vec <- vapply(term2genes_list, function(tg) {
    sum(sample_genes %in% tg)
  }, integer(1), USE.NAMES = FALSE)
  
  # Hypergeometric test
  pvals <- phyper(k_vec - 1, M_vec, N_bg - M_vec, n_target,
                  lower.tail = FALSE)
  
  # BH adjustment
  padj <- p.adjust(pvals, method = "BH")
  
  # Count significant
  iter_sig_counts[i] <- sum(padj < FDR_THRESH, na.rm = TRUE)
  
  if (i %% pb_step == 0) {
    cat(sprintf("    iteration %d/%d (%.0f%%)\n", i, N_ITER, 100*i/N_ITER))
  }
}
t_end <- Sys.time()
elapsed <- difftime(t_end, t_start, units = "secs")
cat(sprintf("  Done in %.1f sec (%.2f sec/iter)\n", elapsed, elapsed/N_ITER))

# Summary statistics
null_median  <- median(iter_sig_counts)
null_q025    <- quantile(iter_sig_counts, 0.025, names = FALSE)
null_q975    <- quantile(iter_sig_counts, 0.975, names = FALSE)
null_mean    <- mean(iter_sig_counts)
null_sd      <- sd(iter_sig_counts)

tier_a_nsig    <- if (is.data.frame(tier_a_res) && nrow(tier_a_res) > 0) nrow(tier_a_res) else 0
near_tc_nsig   <- if (is.data.frame(near_tc_res) && nrow(near_tc_res) > 0) nrow(near_tc_res) else 0
off_only_nsig  <- if (is.data.frame(off_only_res) && nrow(off_only_res) > 0) nrow(off_only_res) else 0

# Z-scores
tier_a_z   <- if (null_sd > 0) (tier_a_nsig - null_mean) / null_sd else NA
near_tc_z  <- if (null_sd > 0) (near_tc_nsig - null_mean) / null_sd else NA

# Empirical p-values
tier_a_p_emp  <- sum(iter_sig_counts >= tier_a_nsig) / N_ITER
near_tc_p_emp <- sum(iter_sig_counts >= near_tc_nsig) / N_ITER

cat(sprintf("\n  Null distribution: median=%d, [%d, %d] (95%% CI)\n",
            null_median, null_q025, null_q975))
cat(sprintf("  Tier A observed: %d (z=%.2f, p_emp=%.4f)\n",
            tier_a_nsig, tier_a_z, tier_a_p_emp))
cat(sprintf("  near-Tier-C observed: %d (z=%.2f, p_emp=%.4f)\n",
            near_tc_nsig, near_tc_z, near_tc_p_emp))

# Save random control results
random_df <- data.frame(
  iteration          = seq_len(N_ITER),
  significant_terms  = iter_sig_counts
)
write_tsv(random_df, file.path(TABLES_DIR, "random_control_results.tsv"))

# ═══════════════════════════════════════════════════════════════════════
# STEP 7: Save enrichment tables
# ═══════════════════════════════════════════════════════════════════════
cat("\n── Step 7: Saving enrichment tables ──\n")

save_enrichment <- function(res, filename) {
  if (is.data.frame(res) && nrow(res) > 0) {
    write_tsv(res, file.path(TABLES_DIR, filename))
    cat(sprintf("  Saved %s (%d rows)\n", filename, nrow(res)))
  } else {
    # Write empty file with header
    writeLines("ID\tDescription\tGeneRatio\tBgRatio\tpvalue\tp.adjust\tqvalue\tgeneID\tCount",
               file.path(TABLES_DIR, filename))
    cat(sprintf("  Saved empty %s\n", filename))
  }
}

save_enrichment(tier_a_res,   "tier_a_enrichment.tsv")
save_enrichment(near_tc_res,  "near_tier_c_enrichment.tsv")
save_enrichment(off_only_res, "off_only_enrichment.tsv")

# ═══════════════════════════════════════════════════════════════════════
# STEP 8: Visualization
# ═══════════════════════════════════════════════════════════════════════
cat("\n── Step 8: Creating visualizations ──\n")

# --- Fig 1: Null distribution histogram ---
cat("  Null distribution histogram...\n")

# Build histogram data
hist_data <- data.frame(sig_terms = iter_sig_counts)

# Determine plot range
x_max <- max(max(iter_sig_counts), tier_a_nsig, near_tc_nsig, off_only_nsig) * 1.15

p_null <- ggplot(hist_data, aes(x = sig_terms)) +
  geom_histogram(bins = 50, fill = "grey70", color = "grey40", alpha = 0.8) +
  # Tier A observed
  geom_vline(xintercept = tier_a_nsig, color = "#E41A1C", linewidth = 1.2,
             linetype = "solid") +
  # near-Tier-C observed
  geom_vline(xintercept = near_tc_nsig, color = "#377EB8", linewidth = 1.2,
             linetype = "dashed") +
  # OFF-only observed
  geom_vline(xintercept = off_only_nsig, color = "#4DAF4A", linewidth = 1.2,
             linetype = "dotted") +
  # Background median
  geom_vline(xintercept = null_median, color = "grey30", linewidth = 0.8,
             linetype = "twodash") +
  annotate("text", x = tier_a_nsig, y = Inf,
           label = sprintf("Tier A: %d", tier_a_nsig),
           color = "#E41A1C", hjust = -0.1, vjust = 2, size = 3.2) +
  annotate("text", x = near_tc_nsig, y = Inf,
           label = sprintf("near-Tier-C: %d", near_tc_nsig),
           color = "#377EB8", hjust = -0.1, vjust = 4, size = 3.2) +
  annotate("text", x = off_only_nsig, y = Inf,
           label = sprintf("OFF-only: %d", off_only_nsig),
           color = "#4DAF4A", hjust = -0.1, vjust = 6, size = 3.2) +
  labs(
    title    = "Random Control: Null Distribution of Significant GO Terms",
    subtitle = sprintf("%d iterations × %d genes (matched to Tier A)",
                       N_ITER, n_target),
    x        = "Number of Significant GO BP Terms (FDR < 0.05)",
    y        = "Frequency",
    caption  = paste0("Vertical lines: observed values | Median: ",
                      null_median)
  ) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(FIGURES_DIR, "random_control_null.pdf"),
       p_null, width = 8, height = 5.5, device = "pdf")
cat("    → random_control_null.pdf\n")

# --- Fig 2: Enrichment comparison barplot ---
cat("  Enrichment comparison barplot...\n")

build_comparison <- function(res, label) {
  if (is.data.frame(res) && nrow(res) > 0) {
    # GeneRatio is stored as "k/n" string; parse it
    gr <- suppressWarnings(as.numeric(sub("/.*", "", res$GeneRatio)) /
                           as.numeric(sub(".*/", "", res$GeneRatio)))
    data.frame(
      set            = label,
      n_sig_terms    = nrow(res),
      median_FDR     = median(res$p.adjust, na.rm = TRUE),
      mean_gene_ratio = mean(gr, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      set             = label,
      n_sig_terms     = NA_real_,
      median_FDR      = NA_real_,
      mean_gene_ratio  = NA_real_,
      stringsAsFactors = FALSE
    )
  }
}

comp_df <- rbind(
  build_comparison(tier_a_res,   "Tier A"),
  build_comparison(near_tc_res,  "near-Tier-C"),
  build_comparison(off_only_res, "OFF-only")
)
comp_df$set <- factor(comp_df$set, levels = c("Tier A", "near-Tier-C", "OFF-only"))

# Reshape for barplot
comp_long <- data.frame(
  set = rep(comp_df$set, 3),
  metric = rep(c("n_sig_terms", "median_FDR", "mean_gene_ratio"), each = 3),
  value = c(comp_df$n_sig_terms, comp_df$median_FDR, comp_df$mean_gene_ratio),
  stringsAsFactors = FALSE
)

p_comp <- ggplot(comp_long, aes(x = set, y = value, fill = set)) +
  geom_col(position = "dodge", width = 0.6) +
  facet_wrap(~ metric, scales = "free_y",
             labeller = labeller(metric = c(
               n_sig_terms      = "Significant Terms",
               median_FDR       = "Median FDR",
               mean_gene_ratio  = "Mean GeneRatio"
             ))) +
  scale_fill_manual(values = c("Tier A" = "#E41A1C",
                               "near-Tier-C" = "#377EB8",
                               "OFF-only" = "#4DAF4A")) +
  labs(title    = "GO Enrichment Comparison Across Gene Sets",
       subtitle = sprintf("Tier A: %d genes, near-Tier-C: %d genes, OFF-only: %d genes",
                          length(tier_a_annot), length(near_tc_annot), length(off_only_annot)),
       x        = NULL,
       y        = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))

ggsave(file.path(FIGURES_DIR, "enrichment_comparison.pdf"),
       p_comp, width = 9, height = 4, device = "pdf")
cat("    → enrichment_comparison.pdf\n")

# --- Fig 3: Top term comparison ---
cat("  Top term comparison dotplot...\n")

# Collect top 10 terms from each set
get_top_terms <- function(res, label, n = 10) {
  if (is.data.frame(res) && nrow(res) > 0) {
    res <- res[order(res$p.adjust), ]
    res <- head(res, n)
    res$set <- label
    res
  } else {
    NULL
  }
}

top_a  <- get_top_terms(tier_a_res,   "Tier A")
top_nt <- get_top_terms(near_tc_res,  "near-Tier-C")
top_of <- get_top_terms(off_only_res, "OFF-only")

top_all <- do.call(rbind, Filter(Negate(is.null), list(top_a, top_nt, top_of)))

if (!is.null(top_all) && nrow(top_all) > 0) {
  top_all$set <- factor(top_all$set, levels = c("Tier A", "near-Tier-C", "OFF-only"))

  # Parse GeneRatio for plotting
  top_all$GeneRatio_num <- suppressWarnings(
    as.numeric(sub("/.*", "", top_all$GeneRatio)) /
    as.numeric(sub(".*/", "", top_all$GeneRatio))
  )

  p_top <- ggplot(top_all, aes(x = GeneRatio_num,
                               y = reorder(Description, GeneRatio_num),
                               size = Count,
                               color = p.adjust)) +
    geom_point(alpha = 0.8) +
    facet_wrap(~ set, scales = "free_y", ncol = 1) +
    scale_color_gradient(low = "red", high = "blue", name = "FDR",
                         trans = "log10") +
    labs(title    = "Top Enriched GO BP Terms by Gene Set",
         subtitle = "Top 10 terms by FDR in each set",
         x        = "GeneRatio",
         y        = NULL,
         size     = "Gene Count") +
    theme(plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"),
          axis.text.y = element_text(size = 8))

  ggsave(file.path(FIGURES_DIR, "top_terms_comparison.pdf"),
         p_top, width = 10, height = max(8, nrow(top_all) * 0.25), device = "pdf")
  cat("    → top_terms_comparison.pdf\n")
}

# ═══════════════════════════════════════════════════════════════════════
# STEP 9: Print final summary
# ═══════════════════════════════════════════════════════════════════════
cat("\n")
cat("═══ ═══ ═══ FINAL SUMMARY ═══ ═══ ═══\n")
cat(sprintf("Tier A:     %d annotated genes, %d significant GO terms (FDR<%.2f)\n",
            length(tier_a_annot), tier_a_nsig, FDR_THRESH))
cat(sprintf("near-Tier-C:%d annotated genes, %d significant GO terms\n",
            length(near_tc_annot), near_tc_nsig))
cat(sprintf("OFF-only:   %d annotated genes, %d significant GO terms\n",
            length(off_only_annot), off_only_nsig))

cat(sprintf("\nRandom control (Tier A size-matched, %d genes):\n", n_target))
cat(sprintf("  Median = %d\n", null_median))
cat(sprintf("  95%% CI = [%d, %d]\n", null_q025, null_q975))
cat(sprintf("  Mean   = %.1f, SD = %.1f\n", null_mean, null_sd))

cat(sprintf("\nTier A vs random control:     z = %.2f, empirical p = %.4f\n",
            tier_a_z, tier_a_p_emp))
cat(sprintf("near-Tier-C vs random control: z = %.2f, empirical p = %.4f\n",
            near_tc_z, near_tc_p_emp))

# Top 5 terms per set
cat("\n── Top 5 GO Terms by FDR ──\n")
print_top5 <- function(res, label) {
  cat(sprintf("\n  %s:\n", label))
  if (is.data.frame(res) && nrow(res) > 0) {
    top <- head(res[order(res$p.adjust), ], 5)
    for (i in seq_len(nrow(top))) {
      cat(sprintf("    %d. %s (FDR=%.2e, Count=%d, GeneRatio=%s)\n",
                  i, substr(top$Description[i], 1, 65),
                  top$p.adjust[i], top$Count[i], top$GeneRatio[i]))
    }
  } else {
    cat("    No significant terms\n")
  }
}
print_top5(tier_a_res,   "Tier A")
print_top5(near_tc_res,  "near-Tier-C")
print_top5(off_only_res, "OFF-only")

cat(sprintf("\nOutput files written to: %s/\n", EXP_DIR))
cat("Done.\n")
