#!/usr/bin/env Rscript
###############################################################################
# OmniQuant-RNA: Bombyx mori Testis vs Ovary — GSEA (Gene Set Enrichment)
#
# Complements the ORA analysis. GSEA uses the full ranked gene list rather
# than threshold-based DEGs, and preserves directionality:
#   rank_score = sign(consensus_logFC) × -log10(min(best_rra_fdr, best_cct_fdr))
#
# Expected biological findings:
#   testis-positive (NES > 0): cilium assembly, dynein, motor proteins, TCA
#   ovary-positive  (NES < 0): ribosome, translation, ribosome biogenesis
###############################################################################

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(GO.db)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
})

source("../theme_nature.R")

# ═════════════════════════════════════════════════════════════════════════════
# 0. Configuration
# ═════════════════════════════════════════════════════════════════════════════

BASE      <- Sys.getenv("PROJECT_ROOT", normalizePath("."))
CONSENSUS <- file.path(BASE,
  "benchmark_results/bombyx_mori/runs",
  "2026-05-24_bombyx_mori_decontam-on",
  "07.consensus_expression/Testis_vs_Ovary",
  "consensus_results.tsv")
GO_MAP <- file.path(BASE,
  "experiments/bombyx_enrichment/data",
  "bombyx_ncbi_to_go_full.tsv")

OUT_DIR <- file.path(BASE, "experiments/bombyx_enrichment/results")
FIG_DIR <- file.path(OUT_DIR, "figures")
TBL_DIR <- file.path(OUT_DIR, "tables")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TBL_DIR, showWarnings = FALSE, recursive = TRUE)

# ═════════════════════════════════════════════════════════════════════════════
# 1. Load GO annotations & build TERM2GENE / TERM2NAME
# ═════════════════════════════════════════════════════════════════════════════

message("=== Loading GO annotations ===")
go_map <- read_tsv(GO_MAP, col_types = cols(.default = col_character()))
message(sprintf("  GO mappings loaded: %d rows, %d unique genes, %d unique GO terms",
                nrow(go_map),
                length(unique(go_map$ncbi_gene_id)),
                length(unique(go_map$GO))))

# TERM2GENE: GO term → NCBI gene ID (unfiltered, used as full universe)
term2gene_all <- go_map %>%
  dplyr::select(GO, ncbi_gene_id) %>%
  dplyr::distinct()
message(sprintf("  TERM2GENE: %d unique mappings", nrow(term2gene_all)))

# Ontology-specific TERM2GENE (filter before each GSEA run)
term2gene_BP <- go_map %>%
  dplyr::filter(ontology == "BP") %>%
  dplyr::select(GO, ncbi_gene_id) %>%
  dplyr::distinct()
term2gene_CC <- go_map %>%
  dplyr::filter(ontology == "CC") %>%
  dplyr::select(GO, ncbi_gene_id) %>%
  dplyr::distinct()
term2gene_MF <- go_map %>%
  dplyr::filter(ontology == "MF") %>%
  dplyr::select(GO, ncbi_gene_id) %>%
  dplyr::distinct()
message(sprintf("  Ontology-specific mappings: BP=%d, CC=%d, MF=%d",
                nrow(term2gene_BP), nrow(term2gene_CC), nrow(term2gene_MF)))

# TERM2NAME: GO term → human-readable description (from GO.db)
go_terms <- unique(go_map$GO)
go_names <- AnnotationDbi::select(GO.db,
  keys    = go_terms,
  columns = c("GOID", "TERM"),
  keytype = "GOID"
)
term2name <- go_names %>%
  dplyr::select(GOID, TERM) %>%
  dplyr::rename(GO = GOID, description = TERM) %>%
  dplyr::distinct()
message(sprintf("  TERM2NAME: %d GO term descriptions", nrow(term2name)))

# ═════════════════════════════════════════════════════════════════════════════
# 2. Load consensus results & build ranked gene list
# ═════════════════════════════════════════════════════════════════════════════

message("\n=== Loading consensus results & building ranked list ===")
consensus <- read_tsv(CONSENSUS, col_types = cols(.default = col_character()))

# Parse gene IDs: "gene:GeneID_692512" → "692512", "gene:LOC101743284" → "101743284"
parse_ncbi <- function(gid) {
  id <- gid
  id <- sub("^gene:GeneID_", "", id)
  id <- sub("^gene:LOC",      "", id)
  id <- sub("^gene:",          "", id)  # fallback for other prefix patterns
  id
}

# Build ranked list
df <- consensus %>%
  dplyr::mutate(
    ncbi_id         = parse_ncbi(gene_id_standard),
    consensus_logFC = as.numeric(consensus_logFC),
    best_rra_fdr    = as.numeric(best_rra_fdr),
    best_cct_fdr    = as.numeric(best_cct_fdr),

    # Take the better (smaller) FDR between RRA and CCT consensus
    # pmin(..., na.rm = TRUE) handles NA: returns the non-NA value
    best_fdr = pmin(best_rra_fdr, best_cct_fdr, na.rm = TRUE),

    # Safety: cap at 1e-300 to avoid -log10(0) = Inf
    # (consensus pipeline already clips at 1e-16, but this guards edge cases)
    best_fdr = pmax(best_fdr, 1e-300),

    dir_sign   = sign(consensus_logFC),
    rank_score = dir_sign * (-log10(best_fdr))
  ) %>%
  # Filter: remove genes without valid logFC or FDR
  dplyr::filter(!is.na(rank_score), is.finite(rank_score)) %>%
  dplyr::filter(!is.na(ncbi_id), ncbi_id != "")

message(sprintf("  Genes with valid rank_score: %d", nrow(df)))
message(sprintf("    Positive (testis-high, rank > 0): %d", sum(df$rank_score > 0)))
message(sprintf("    Negative (ovary-high,  rank < 0): %d", sum(df$rank_score < 0)))
message(sprintf("    Rank score range: [%.2f, %.2f]",
                min(df$rank_score), max(df$rank_score)))

# Build named sorted vector for GSEA
gene_list <- df$rank_score
names(gene_list) <- df$ncbi_id
gene_list <- sort(gene_list, decreasing = TRUE)

# Remove duplicate gene names — keep first (highest rank_score)
gene_list <- gene_list[!duplicated(names(gene_list))]
message(sprintf("  Final ranked list (unique genes): %d", length(gene_list)))

# ═════════════════════════════════════════════════════════════════════════════
# 3. Run GSEA for BP, CC, MF (ontology-specific TERM2GENE)
# ═════════════════════════════════════════════════════════════════════════════

run_gsea <- function(gene_list, term2gene, term2name, ontology_label) {
  message(sprintf("\n--- Running GSEA (%s) ---", ontology_label))
  set.seed(42)
  res <- GSEA(
    geneList      = gene_list,
    TERM2GENE     = term2gene,
    TERM2NAME     = term2name,
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    minGSSize     = 5,
    maxGSSize     = 500,
    seed          = TRUE,
    verbose       = FALSE
  )
  n_sig <- if (nrow(res) > 0) sum(res@result$p.adjust < 0.05) else 0
  n_pos <- if (nrow(res) > 0) sum(res@result$p.adjust < 0.05 & res@result$NES > 0) else 0
  n_neg <- if (nrow(res) > 0) sum(res@result$p.adjust < 0.05 & res@result$NES < 0) else 0
  message(sprintf("  Total terms: %d, Significant (FDR<0.05): %d", nrow(res), n_sig))
  message(sprintf("    Testis-positive (NES>0): %d", n_pos))
  message(sprintf("    Ovary-positive  (NES<0): %d", n_neg))
  res
}

gsea_bp <- run_gsea(gene_list, term2gene_BP, term2name, "BP")
gsea_cc <- run_gsea(gene_list, term2gene_CC, term2name, "CC")
gsea_mf <- run_gsea(gene_list, term2gene_MF, term2name, "MF")

# ═════════════════════════════════════════════════════════════════════════════
# 4. Save results as TSV
# ═════════════════════════════════════════════════════════════════════════════

message("\n=== Saving GSEA results ===")
gsea_results <- list(BP = gsea_bp, CC = gsea_cc, MF = gsea_mf)

for (ont in names(gsea_results)) {
  res <- gsea_results[[ont]]
  if (nrow(res) > 0) {
    out_path <- file.path(TBL_DIR, paste0("gsea_results_", ont, ".tsv"))
    write_tsv(as.data.frame(res), out_path)
    message(sprintf("  Saved %s (%d rows)", out_path, nrow(res)))
  }
}

# Also save directional subsets
if (nrow(gsea_bp) > 0) {
  bp_pos <- gsea_bp@result %>%
    dplyr::filter(NES > 0, p.adjust < 0.05) %>%
    dplyr::arrange(desc(NES))
  bp_neg <- gsea_bp@result %>%
    dplyr::filter(NES < 0, p.adjust < 0.05) %>%
    dplyr::arrange(NES)

  write_tsv(bp_pos, file.path(TBL_DIR, "gsea_bp_testis_positive.tsv"))
  write_tsv(bp_neg, file.path(TBL_DIR, "gsea_bp_ovary_positive.tsv"))
  message(sprintf("  Testis-positive (FDR<0.05): %d terms", nrow(bp_pos)))
  message(sprintf("  Ovary-positive  (FDR<0.05): %d terms", nrow(bp_neg)))
}

# ═════════════════════════════════════════════════════════════════════════════
# 5. GSEA curve plots
# ═════════════════════════════════════════════════════════════════════════════

message("\n=== Generating GSEA curve plots ===")

# Top testis-positive BP terms (NES > 0, FDR < 0.05)
top_testis <- tryCatch({
  gsea_bp@result %>%
    dplyr::filter(NES > 0, p.adjust < 0.05) %>%
    head(4)
}, error = function(e) NULL)

# Top ovary-positive BP terms (NES < 0, FDR < 0.05)
top_ovary <- tryCatch({
  gsea_bp@result %>%
    dplyr::filter(NES < 0, p.adjust < 0.05) %>%
    head(4)
}, error = function(e) NULL)

# Testis-positive curves (NES > 0)
if (!is.null(top_testis) && nrow(top_testis) > 0) {
  pdf(file.path(FIG_DIR, "gsea_testis_positive_curves.pdf"),
      width = 12, height = 10)
  for (i in seq_len(nrow(top_testis))) {
    term_id <- rownames(top_testis)[i]
    p <- gseaplot2(gsea_bp,
                   geneSetID    = term_id,
                   title        = top_testis$Description[i],
                   pvalue_table = FALSE)
    print(p)
  }
  dev.off()
  message(sprintf("  Saved gsea_testis_positive_curves.pdf (%d terms)", nrow(top_testis)))
} else {
  warning("  No testis-positive terms with FDR<0.05 for curves")
}

# Ovary-positive curves (NES < 0)
if (!is.null(top_ovary) && nrow(top_ovary) > 0) {
  pdf(file.path(FIG_DIR, "gsea_ovary_positive_curves.pdf"),
      width = 12, height = 10)
  for (i in seq_len(nrow(top_ovary))) {
    term_id <- rownames(top_ovary)[i]
    p <- gseaplot2(gsea_bp,
                   geneSetID    = term_id,
                   title        = top_ovary$Description[i],
                   pvalue_table = FALSE)
    print(p)
  }
  dev.off()
  message(sprintf("  Saved gsea_ovary_positive_curves.pdf (%d terms)", nrow(top_ovary)))
} else {
  warning("  No ovary-positive terms with FDR<0.05 for curves")
}

# Ridgeplot: top 20 BP terms by adjusted p-value
if (nrow(gsea_bp) > 0) {
  pdf(file.path(FIG_DIR, "gsea_bp_ridgeplot.pdf"), width = 12, height = 16)
  print(
    ridgeplot(gsea_bp, showCategory = 20, fill = "p.adjust") +
      labs(title = "GSEA: GO Biological Process — Testis vs Ovary")
  )
  dev.off()
  message("  Saved gsea_bp_ridgeplot.pdf")
}

# ═════════════════════════════════════════════════════════════════════════════
# 6. Summary
# ═════════════════════════════════════════════════════════════════════════════

message("\n=== GSEA SUMMARY ===")

for (ont in names(gsea_results)) {
  res <- gsea_results[[ont]]
  if (nrow(res) == 0) next
  sig <- dplyr::filter(res@result, p.adjust < 0.05)
  n_pos <- sum(sig$NES > 0)
  n_neg <- sum(sig$NES < 0)
  message(sprintf("\n--- %s ---", ont))
  message(sprintf("  Significant terms (FDR<0.05): %d total  (%d testis-positive, %d ovary-positive)",
                  nrow(sig), n_pos, n_neg))
}

# Top 5 testis-positive BP terms
message("\n=== Top 5 Testis-Positive BP Terms (NES > 0, FDR < 0.05) ===")
if (nrow(gsea_bp) > 0) {
  bp_pos5 <- gsea_bp@result %>%
    dplyr::filter(NES > 0, p.adjust < 0.05) %>%
    head(5)
  for (i in seq_len(nrow(bp_pos5))) {
    message(sprintf("  %d. %s  (NES=%.2f, FDR=%.2e, core=%d/%d)",
                    i,
                    bp_pos5$Description[i],
                    bp_pos5$NES[i],
                    bp_pos5$p.adjust[i],
                    length(strsplit(bp_pos5$core_enrichment[i], "/")[[1]]),
                    bp_pos5$setSize[i]))
  }
} else {
  message("  (none)")
}

# Top 5 ovary-positive BP terms
message("\n=== Top 5 Ovary-Positive BP Terms (NES < 0, FDR < 0.05) ===")
if (nrow(gsea_bp) > 0) {
  bp_neg5 <- gsea_bp@result %>%
    dplyr::filter(NES < 0, p.adjust < 0.05) %>%
    head(5)
  for (i in seq_len(nrow(bp_neg5))) {
    message(sprintf("  %d. %s  (NES=%.2f, FDR=%.2e, core=%d/%d)",
                    i,
                    bp_neg5$Description[i],
                    bp_neg5$NES[i],
                    bp_neg5$p.adjust[i],
                    length(strsplit(bp_neg5$core_enrichment[i], "/")[[1]]),
                    bp_neg5$setSize[i]))
  }
} else {
  message("  (none)")
}

# ═════════════════════════════════════════════════════════════════════════════
# 7. Save full results object
# ═════════════════════════════════════════════════════════════════════════════

saveRDS(
  list(BP = gsea_bp, CC = gsea_cc, MF = gsea_mf, ranked_list = gene_list),
  file.path(OUT_DIR, "gsea_results.rds")
)
message("\n=== GSEA COMPLETE ===")
message(sprintf("  RDS:  %s", file.path(OUT_DIR, "gsea_results.rds")))
message(sprintf("  TSVs: %s/", TBL_DIR))
message(sprintf("  PDFs: %s/", FIG_DIR))
