#!/usr/bin/env Rscript
###############################################################################
# OmniQuant-RNA: Bombyx mori Testis vs Ovary — GO/KEGG Enrichment Analysis
#
# Compares consensus Tier A (up/down) vs Tier C vs background for 
# functional enrichment analysis. Demonstrates that Tier A genes are
# enriched in tissue-specific biological processes while Tier C genes
# are enriched in housekeeping/noise pathways.
###############################################################################

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
})

# ── Configuration ────────────────────────────────────────────────────────────

DATA_DIR   <- "experiments/bombyx_enrichment/data"
GENESET_DIR <- "experiments/bombyx_enrichment/genesets"
OUT_DIR    <- "experiments/bombyx_enrichment/results"
FIG_DIR    <- file.path(OUT_DIR, "figures")
TBL_DIR    <- file.path(OUT_DIR, "tables")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TBL_DIR, recursive = TRUE, showWarnings = FALSE)

# ── 1. Load Annotations ──────────────────────────────────────────────────────

message("Loading annotation data...")

# GO annotations (UniProt-GOA via NCBI GeneID)
go_map <- read_tsv(
  file.path(DATA_DIR, "bombyx_ncbi_to_go_full.tsv"),
  col_types = cols(.default = col_character())
)

# Build TERM2GENE: GO term → NCBI Gene ID
term2gene <- go_map %>% dplyr::select(GO, ncbi_gene_id) %>% dplyr::distinct()
message(sprintf("  TERM2GENE: %d unique GO-term-to-gene mappings", nrow(term2gene)))

# GO term descriptions (from QuickGO or GO.db)
# We need TERM2NAME for the enricher
# Since we don't have GO descriptions in the TSV, we'll use GO.db
if (requireNamespace("GO.db", quietly = TRUE)) {
  library(GO.db)
  go_terms <- unique(go_map$GO)
  go_names <- AnnotationDbi::select(GO.db, keys = go_terms, columns = c("GOID", "TERM"), keytype = "GOID")
  term2name <- go_names %>% dplyr::select(GOID, TERM) %>% 
    dplyr::rename(GO = GOID, description = TERM)
} else {
  # Fallback: use GO ID as name
  term2name <- data.frame(
    GO = unique(go_map$GO),
    description = unique(go_map$GO),
    stringsAsFactors = FALSE
  )
}
message(sprintf("  TERM2NAME: %d GO term descriptions", nrow(term2name)))

# ── 2. Load Gene Sets ────────────────────────────────────────────────────────

message("Loading gene sets...")

load_gene_set <- function(name) {
  path <- file.path(GENESET_DIR, paste0(name, ".txt"))
  if (!file.exists(path)) stop(sprintf("File not found: %s", path))
  genes <- readLines(path)
  genes <- genes[genes != ""]
  message(sprintf("  %s: %d genes", name, length(genes)))
  genes
}

tier_a_up     <- load_gene_set("tier_a_up")
tier_a_down   <- load_gene_set("tier_a_down")
tier_c        <- load_gene_set("tier_c")
background    <- load_gene_set("unclassified_no_signal")
all_bg        <- load_gene_set("all_with_go")

# ── 3. GO Enrichment ─────────────────────────────────────────────────────────

message("Running GO enrichment (BP)...")

run_go_enrichment <- function(gene_list, bg_list, ont = "BP", name = "gene_set") {
  # Filter TERM2GENE to the specific ontology
  go_terms_bp <- go_map %>% 
    dplyr::filter(ontology == ont) %>% 
    dplyr::select(GO, ncbi_gene_id) %>% 
    dplyr::distinct()
  
  if (nrow(go_terms_bp) == 0) {
    warning(sprintf("No %s terms found", ont))
    return(NULL)
  }
  
  res <- enricher(
    gene     = gene_list,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    minGSSize = 5,
    maxGSSize = 500,
    qvalueCutoff = 0.2,
    universe  = bg_list,
    TERM2GENE = go_terms_bp,
    TERM2NAME = term2name
  )
  
  if (!is.null(res) && nrow(res) > 0) {
    res@result$GeneSet <- name
    res@result$Ontology <- ont
  }
  res
}

# Run enrichment for each comparison group
go_results <- list()

# Tier A up vs all background
go_results$tier_a_up_BP <- run_go_enrichment(tier_a_up, all_bg, "BP", "Tier_A_up")
go_results$tier_a_up_CC <- run_go_enrichment(tier_a_up, all_bg, "CC", "Tier_A_up")
go_results$tier_a_up_MF <- run_go_enrichment(tier_a_up, all_bg, "MF", "Tier_A_up")

# Tier A down vs all background
go_results$tier_a_down_BP <- run_go_enrichment(tier_a_down, all_bg, "BP", "Tier_A_down")
go_results$tier_a_down_CC <- run_go_enrichment(tier_a_down, all_bg, "CC", "Tier_A_down")
go_results$tier_a_down_MF <- run_go_enrichment(tier_a_down, all_bg, "MF", "Tier_A_down")

# Tier C vs all background (control)
go_results$tier_c_BP <- run_go_enrichment(tier_c, all_bg, "BP", "Tier_C")
go_results$tier_c_CC <- run_go_enrichment(tier_c, all_bg, "CC", "Tier_C")
go_results$tier_c_MF <- run_go_enrichment(tier_c, all_bg, "MF", "Tier_C")

# Summarize
for (name in names(go_results)) {
  res <- go_results[[name]]
  if (!is.null(res) && nrow(res) > 0) {
    n_sig <- sum(res@result$p.adjust < 0.05)
    message(sprintf("  %s: %d terms, %d significant (FDR<0.05)", name, nrow(res), n_sig))
  }
}

# ── 4. Save GO Results ──────────────────────────────────────────────────────

message("Saving GO enrichment results...")

for (name in names(go_results)) {
  res <- go_results[[name]]
  if (!is.null(res) && nrow(res) > 0) {
    out_path <- file.path(TBL_DIR, paste0("go_enrichment_", name, ".tsv"))
    write_tsv(as.data.frame(res), out_path)
    message(sprintf("  Saved %s", out_path))
  }
}

# ── 5. GO Visualization ──────────────────────────────────────────────────────

message("Generating GO enrichment plots...")

plot_go_dotplot <- function(results_list, title, outfile, top_n = 15) {
  # Combine multiple enrichment results into one dotplot
  valid <- results_list[!sapply(results_list, is.null)]
  valid <- valid[sapply(valid, function(x) nrow(x) > 0)]
  
  if (length(valid) == 0) {
    message(sprintf("  Skipping dotplot '%s': no results", title))
    return(invisible())
  }
  
  # Combine and keep top N per set
  combined <- bind_rows(lapply(names(valid), function(nm) {
    res <- valid[[nm]]
    df <- as.data.frame(res)
    if (nrow(df) > top_n) df <- df[1:top_n, ]
    df
  }))
  
  p <- dotplot(valid[[1]], showCategory = top_n, title = title) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(size = 13, face = "bold"),
      axis.text.y = element_text(size = 9)
    )
  
  ggsave(outfile, p, width = 10, height = 8, dpi = 300)
  message(sprintf("  Saved dotplot: %s", outfile))
}

# Generate key plots
if (!is.null(go_results$tier_a_up_BP) && nrow(go_results$tier_a_up_BP) > 0) {
  pdf(file.path(FIG_DIR, "go_tier_a_up_BP_dotplot.pdf"), width = 10, height = 8)
  print(dotplot(go_results$tier_a_up_BP, showCategory = 20, 
                title = "Tier A Up (Testis-enriched): GO Biological Process"))
  dev.off()
  message("  Saved Tier A up BP dotplot")
}

if (!is.null(go_results$tier_a_down_BP) && nrow(go_results$tier_a_down_BP) > 0) {
  pdf(file.path(FIG_DIR, "go_tier_a_down_BP_dotplot.pdf"), width = 10, height = 8)
  print(dotplot(go_results$tier_a_down_BP, showCategory = 20,
                title = "Tier A Down (Ovary-enriched): GO Biological Process"))
  dev.off()
  message("  Saved Tier A down BP dotplot")
}

if (!is.null(go_results$tier_c_BP) && nrow(go_results$tier_c_BP) > 0) {
  pdf(file.path(FIG_DIR, "go_tier_c_BP_dotplot.pdf"), width = 10, height = 8)
  print(dotplot(go_results$tier_c_BP, showCategory = 20,
                title = "Tier C: GO Biological Process"))
  dev.off()
  message("  Saved Tier C BP dotplot")
}

# Side-by-side comparison: Tier A up vs Tier C
if (!is.null(go_results$tier_a_up_BP) && nrow(go_results$tier_a_up_BP) > 0 &&
    !is.null(go_results$tier_c_BP) && nrow(go_results$tier_c_BP) > 0) {
  
  # Get top terms from each
  a_up_terms <- head(go_results$tier_a_up_BP@result$Description, 20)
  c_terms    <- head(go_results$tier_c_BP@result$Description, 20)
  
  # Overlap analysis
  overlap <- intersect(a_up_terms, c_terms)
  a_unique <- setdiff(a_up_terms, c_terms)
  c_unique <- setdiff(c_terms, a_up_terms)
  
  message(sprintf("  Tier A up vs Tier C overlap: %d shared terms", length(overlap)))
  message(sprintf("  Tier A up unique: %d terms", length(a_unique)))
  message(sprintf("  Tier C unique: %d terms", length(c_unique)))
  
  # Save overlap
  writeLines(c(
    sprintf("# Tier A up unique terms (%d)", length(a_unique)),
    a_unique,
    "",
    sprintf("# Tier C unique terms (%d)", length(c_unique)),
    c_unique,
    "",
    sprintf("# Shared terms (%d)", length(overlap)),
    overlap
  ), file.path(TBL_DIR, "go_term_comparison_tierA_vs_tierC.txt"))
}

# ── 6. KEGG Pathway Enrichment ───────────────────────────────────────────────

message("Running KEGG pathway enrichment...")

# KEGG for Bombyx mori uses organism code "bmor"
# enrichKEGG requires NCBI Gene IDs as input
run_kegg <- function(gene_list, bg_list, name = "gene_set") {
  tryCatch({
    res <- enrichKEGG(
      gene     = gene_list,
      organism = "bmor",
      keyType  = "ncbi-geneid",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      minGSSize = 5,
      maxGSSize = 500,
      qvalueCutoff = 0.2,
      universe  = bg_list
    )
    if (!is.null(res) && nrow(res) > 0) {
      message(sprintf("  %s (KEGG): %d pathways, %d significant", 
                       name, nrow(res), sum(res@result$p.adjust < 0.05)))
    }
    res
  }, error = function(e) {
    message(sprintf("  %s (KEGG): ERROR - %s", name, e$message))
    NULL
  })
}

kegg_results <- list()
kegg_results$tier_a_up   <- run_kegg(tier_a_up, all_bg, "Tier_A_up")
kegg_results$tier_a_down <- run_kegg(tier_a_down, all_bg, "Tier_A_down")
kegg_results$tier_c      <- run_kegg(tier_c, all_bg, "Tier_C")

# Save KEGG results
for (name in names(kegg_results)) {
  res <- kegg_results[[name]]
  if (!is.null(res) && nrow(res) > 0) {
    write_tsv(as.data.frame(res), file.path(TBL_DIR, paste0("kegg_enrichment_", name, ".tsv")))
    
    if (nrow(res) >= 3) {
      pdf(file.path(FIG_DIR, paste0("kegg_", name, "_dotplot.pdf")), width = 10, height = 7)
      print(dotplot(res, showCategory = 15, title = paste(name, ": KEGG Pathways")))
      dev.off()
    }
  }
}

# ── 7. Summary Report ────────────────────────────────────────────────────────

message("\n=== ENRICHMENT ANALYSIS COMPLETE ===")
message(sprintf("GO BP results: %d gene sets analyzed", 
                sum(!sapply(go_results[grepl("BP$", names(go_results))], is.null))))
message(sprintf("KEGG results: %d gene sets analyzed", 
                sum(!sapply(kegg_results, is.null))))
message(sprintf("Figures saved to: %s", FIG_DIR))
message(sprintf("Tables saved to: %s", TBL_DIR))

# ── 8. Key Numbers for Paper ─────────────────────────────────────────────────

message("\n=== KEY NUMBERS FOR PAPER ===")

# Tier A up specific pathways
if (!is.null(go_results$tier_a_up_BP) && nrow(go_results$tier_a_up_BP) > 0) {
  top5 <- head(go_results$tier_a_up_BP@result, 5)
  message("\nTop 5 GO BP terms for Tier A up (testis-enriched):")
  for (i in 1:nrow(top5)) {
    message(sprintf("  %d. %s (FDR=%.2e, %d/%d genes)", 
                    i, top5$Description[i], top5$p.adjust[i], 
                    top5$Count[i], length(tier_a_up)))
  }
}

# Tier C enriched pathways (expected: housekeeping)
if (!is.null(go_results$tier_c_BP) && nrow(go_results$tier_c_BP) > 0) {
  top5 <- head(go_results$tier_c_BP@result, 5)
  message("\nTop 5 GO BP terms for Tier C (low-confidence control):")
  for (i in 1:nrow(top5)) {
    message(sprintf("  %d. %s (FDR=%.2e, %d/%d genes)", 
                    i, top5$Description[i], top5$p.adjust[i], 
                    top5$Count[i], length(tier_c)))
  }
}

message("\nDone.")
