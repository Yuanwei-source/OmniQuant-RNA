#!/usr/bin/env Rscript
###############################################################################
# Figure 6: Functional coherence of consensus tiers
# Bombyx testis–ovary transcriptomes — 6-panel publication-quality figure
#
# Panels:
#   A = Validation design schematic (new)
#   B = Tier A up GO/KEGG dotplot (from PDF)
#   C = Tier A down GO/KEGG dotplot (from PDF)
#   D = GSEA enrichment curves (from PDFs)
#   E = Coherence score comparison barplot
#   F = Not-in-FC functional enrichment dotplot
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(pdftools)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(clusterProfiler)
  library(grid)
  library(gridExtra)
})

source("paper_analysis/theme_nature.R")

# ── Paths ─────────────────────────────────────────────────────────────────────
PROJECT  <- Sys.getenv("PROJECT_ROOT", normalizePath("."))
BOMBYX_FIG  <- file.path(PROJECT, "experiments/bombyx_enrichment/results/figures")
BOMBYX_TBL  <- file.path(PROJECT, "experiments/bombyx_enrichment/results/tables")
BOMBYX_GENE <- file.path(PROJECT, "experiments/bombyx_enrichment/genesets")
BOMBYX_DATA <- file.path(PROJECT, "experiments/bombyx_enrichment/data")
DROS_FIG    <- file.path(PROJECT, "experiments/drosophila_enrichment/results/figures")
DROS_TBL    <- file.path(PROJECT, "experiments/drosophila_enrichment/results/tables")
CONSENSUS   <- file.path(PROJECT,
  "benchmark_results/bombyx_mori/runs/2026-05-24_bombyx_mori_decontam-on",
  "07.consensus_expression/Testis_vs_Ovary/consensus_results.tsv")
BENCH_FIG   <- file.path(PROJECT, "benchmark_results/figures")
EXP_FIG     <- file.path(PROJECT, "experiments/bombyx_enrichment/results/figures")

dir.create(BENCH_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(EXP_FIG,  recursive = TRUE, showWarnings = FALSE)

# ── Global theme ──────────────────────────────────────────────────────────────
theme_pub <- theme_nature(base_size = 6.5, base_family = "Arial") +
  theme(
    axis.line         = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks        = element_line(linewidth = 0.3, colour = "black"),
    panel.grid        = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.15, colour = "grey90"),
    legend.position   = "right",
    legend.title      = element_text(size = 7),
    legend.text       = element_text(size = 6.5),
    plot.title        = element_text(size = 9, face = "bold", hjust = 0.5),
    plot.subtitle     = element_text(size = 7, hjust = 0.5, colour = "grey40"),
    strip.text        = element_text(size = 7, face = "bold"),
    axis.title        = element_text(size = 7.5),
    axis.text         = element_text(size = 7),
    plot.margin       = margin(8, 8, 6, 6)
  )

# ── Helper: add panel label ───────────────────────────────────────────────────
add_label <- function(p, label) {
  p + labs(tag = label) +
    theme(
      plot.tag          = element_text(size = 12, face = "bold", family = ""),
      plot.tag.position = "topleft"
    )
}

# ── Helper: render PDF page as raster for ggplot ──────────────────────────────
pdf_as_gg <- function(pdf_path, page = 1, dpi = 300) {
  tmp <- tempfile(fileext = ".png")
  bitmap <- pdftools::pdf_render_page(pdf_path, page = page, dpi = dpi)
  png::writePNG(bitmap, tmp)
  img <- png::readPNG(tmp)
  unlink(tmp)
  p <- ggdraw() + draw_image(img)
  return(p)
}

# ═══════════════════════════════════════════════════════════════════════════════
# PANEL A: Validation Design Schematic
# ═══════════════════════════════════════════════════════════════════════════════

make_panel_a <- function() {
  # Data frame for box positions (x, y, width, height, label, description, count, color)
  boxes <- data.frame(
    x     = c(0.05, 0.38, 0.72, 0.72, 0.72, 0.72),
    y     = c(0.55, 0.55, 0.80, 0.58, 0.36, 0.14),
    w     = c(0.28, 0.28, 0.23, 0.23, 0.23, 0.23),
    h     = c(0.20, 0.20, 0.17, 0.17, 0.17, 0.17),
    label = c("Consensus\nPipeline", "Tier\nClassification", "Tier A up\n(testis-enriched)",
              "Tier A down\n(ovary-enriched)", "Tier C\n(low-confidence)",
              "OFF-only\n(decontam-sensitive)"),
    sublabel = c("4 quantifiers\nDESeq2 → RRA + CCT", "Evidence-based\ntier assignment", 
                 "n = 2,693\ntop: cilium assembly", "n = 1,546\ntop: translation",
                 "n = 507\nhousekeeping terms", "n = 384\n0 significant terms"),
    fill  = c("#333333", "#666666", "#E74C3C", "#3498DB", "#95A5A6", "#BDC3C7"),
    group = c("input", "input", "group", "group", "group", "group")
  )
  
  # Arrows between boxes
  arrows <- data.frame(
    x    = c(0.33, 0.66),
    y    = c(0.65, 0.65),
    xend = c(0.38, 0.72),
    yend = c(0.65, 0.65)
  )
  
  # Vertical bracket for the 5 groups
  bracket <- data.frame(
    x = c(0.94, 0.97, 0.97),
    y = c(0.14, 0.14, 0.88),
    group = c(1,1,1)
  )
  
  p <- ggplot() +
    # Draw boxes
    geom_rect(data = boxes,
              aes(xmin = x, xmax = x + w,
                  ymin = y - h/2, ymax = y + h/2,
                  fill = fill),
              colour = "white", linewidth = 1.2, alpha = 0.92) +
    # Box labels
    geom_text(data = boxes,
              aes(x = x + w/2, y = y + 0.035,
                  label = label, colour = ifelse(group == "input", "white", "white")),
              size = 3.0, fontface = "bold", lineheight = 0.9, family = "") +
    geom_text(data = boxes,
              aes(x = x + w/2, y = y - 0.045,
                  label = sublabel, colour = ifelse(group == "input", "#CCCCCC", "#F0F0F0")),
              size = 2.3, lineheight = 0.9, family = "") +
    # Arrows
    geom_segment(data = arrows,
                 aes(x = x, y = y, xend = xend, yend = yend),
                 arrow = arrow(length = unit(3, "mm"), type = "closed"),
                 linewidth = 1.0, colour = "#444444") +
    # Section labels
    annotate("text", x = 0.19, y = 0.87, label = "Quantitative\nIntegration",
             size = 2.8, fontface = "italic", colour = "#666666", family = "") +
    annotate("text", x = 0.52, y = 0.87, label = "Gene\nAssignment",
             size = 2.8, fontface = "italic", colour = "#666666", family = "") +
    annotate("text", x = 0.835, y = 0.87, label = "Biological\nValidation",
             size = 2.8, fontface = "italic", colour = "#666666", family = "") +
    # Header
    annotate("text", x = 0.5, y = 0.97, 
             label = "Validation Design: 5 Comparison Groups",
             size = 3.5, fontface = "bold", family = "", colour = "#222222") +
    scale_fill_identity() +
    scale_colour_identity() +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void() +
    theme(plot.margin = margin(5, 10, 5, 10))
  
  add_label(p, "A")
}

# ═══════════════════════════════════════════════════════════════════════════════
# PANEL B: Tier A up GO/KEGG Dotplot (extract from PDF)
# ═══════════════════════════════════════════════════════════════════════════════

make_panel_b <- function() {
  pdf_path <- file.path(BOMBYX_FIG, "go_tier_a_up_BP_dotplot.pdf")
  p <- pdf_as_gg(pdf_path)
  
  # Wrap in ggplot with title
  p_title <- ggplot() + 
    annotate("text", x = 0.5, y = 0.5, 
             label = "Tier A up (testis-enriched): GO Biological Process",
             size = 3.2, fontface = "bold", family = "") +
    theme_void()
  
  combined <- p_title / p + plot_layout(heights = c(0.08, 1))
  add_label(combined, "B")
}

# ═══════════════════════════════════════════════════════════════════════════════
# PANEL C: Tier A down GO/KEGG Dotplot (extract from PDF)
# ═══════════════════════════════════════════════════════════════════════════════

make_panel_c <- function() {
  pdf_path <- file.path(BOMBYX_FIG, "go_tier_a_down_BP_dotplot.pdf")
  p <- pdf_as_gg(pdf_path)
  
  p_title <- ggplot() + 
    annotate("text", x = 0.5, y = 0.5, 
             label = "Tier A down (ovary-enriched): GO Biological Process",
             size = 3.2, fontface = "bold", family = "") +
    theme_void()
  
  combined <- p_title / p + plot_layout(heights = c(0.08, 1))
  add_label(combined, "C")
}

# ═══════════════════════════════════════════════════════════════════════════════
# PANEL D: GSEA Enrichment Curves
# ═══════════════════════════════════════════════════════════════════════════════

make_panel_d <- function() {
  # Extract both GSEA curve PDFs
  p_testis <- pdf_as_gg(file.path(BOMBYX_FIG, "gsea_testis_positive_curves.pdf"), page = 1)
  p_ovary  <- pdf_as_gg(file.path(BOMBYX_FIG, "gsea_ovary_positive_curves.pdf"), page = 1)
  
  # Create title strips
  p1_label <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = "GSEA: Testis-positive (cilium assembly, NES=+2.60)",
             size = 3.2, fontface = "bold", family = "", colour = "#E74C3C") +
    theme_void()
  
  p2_label <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = "GSEA: Ovary-positive (translation, NES=−4.43)",
             size = 3.2, fontface = "bold", family = "", colour = "#3498DB") +
    theme_void()
  
  # Arrange side-by-side
  top_row  <- p1_label + p2_label + plot_layout(widths = c(1, 1))
  bot_row  <- p_testis + p_ovary + plot_layout(widths = c(1, 1))
  combined <- top_row / bot_row + plot_layout(heights = c(0.08, 0.92))
  
  add_label(combined, "D")
}

# ═══════════════════════════════════════════════════════════════════════════════
# PANEL E: Coherence Score Comparison
# ═══════════════════════════════════════════════════════════════════════════════

make_panel_e <- function() {
  # Data from verified enrichment counts
  df <- data.frame(
    Group = c("Tier A\n(Bombyx)", "Tier A\n(Drosophila)", 
              "near-Tier-C\n(Drosophila)", "OFF-only\n(Drosophila)",
              "Random\n(median)"),
    Terms = c(16, 57, 2, 0, 0),
    Color = c("#E74C3C", "#C0392B", "#E67E22", "#BDC3C7", "#95A5A6"),
    Category = c("Tier A", "Tier A", "Control", "Control", "Control")
  )
  df$Group <- factor(df$Group, levels = rev(df$Group))
  
  p <- ggplot(df, aes(x = Group, y = Terms, fill = Color)) +
    geom_col(width = 0.65, colour = NA) +
    geom_text(aes(label = sprintf("%d", Terms), y = Terms + max(Terms)*0.04),
              size = 3.5, fontface = "bold", family = "", colour = "#333333") +
    scale_fill_identity() +
    labs(
      title = "Enrichment coherence across tiers & species",
      subtitle = "Number of significant GO BP terms (FDR<0.05)",
      x = NULL, y = "Significant GO BP terms"
    ) +
    coord_flip() +
    theme_pub +
    theme(
      axis.text.y = element_text(size = 8, lineheight = 0.9),
      panel.grid.major.x = element_line(linewidth = 0.15, colour = "grey85")
    )
  
  add_label(p, "E")
}

# ═══════════════════════════════════════════════════════════════════════════════
# PANEL F: Not-in-FC Functional Enrichment
# ═══════════════════════════════════════════════════════════════════════════════

make_panel_f <- function() {
  message("  Panel F: Running enrichment for not-in-FC genes...")
  
  # Load TERM2GENE
  go_map <- read_tsv(
    file.path(BOMBYX_DATA, "bombyx_ncbi_to_go_full.tsv"),
    col_types = cols(.default = col_character()),
    show_col_types = FALSE
  )
  
  go_terms_bp <- go_map %>%
    dplyr::filter(ontology == "BP") %>%
    dplyr::select(GO, ncbi_gene_id) %>%
    dplyr::distinct()
  
  # TERM2NAME: use GO IDs as names (simple approach)
  term2name <- data.frame(
    GO = unique(go_map$GO),
    description = unique(go_map$GO),
    stringsAsFactors = FALSE
  )
  message(sprintf("    TERM2GENE: %d rows for BP", nrow(go_terms_bp)))
  
  # Load gene sets
  tier_a_not_fc <- readLines(file.path(BOMBYX_GENE, "tier_a_not_fc.txt"))
  tier_a_not_fc <- tier_a_not_fc[tier_a_not_fc != ""]
  background    <- readLines(file.path(BOMBYX_GENE, "all_with_go.txt"))
  background    <- background[background != ""]
  
  message(sprintf("    Query genes: %d, Background: %d", length(tier_a_not_fc), length(background)))
  
  # Run enrichment
  res <- enricher(
    gene     = tier_a_not_fc,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    minGSSize = 5,
    maxGSSize = 500,
    qvalueCutoff = 0.2,
    universe  = background,
    TERM2GENE = go_terms_bp,
    TERM2NAME = term2name
  )
  
  if (is.null(res) || nrow(res) == 0) {
    warning("Panel F: No enrichment results. Creating placeholder.")
    p <- ggplot() + 
      annotate("text", x = 0.5, y = 0.5, label = "No significant enrichment\nfor not-in-FC genes",
               size = 4, family = "") +
      theme_void()
  } else {
    n_sig <- sum(res@result$p.adjust < 0.05)
    message(sprintf("    Significant terms: %d / %d", n_sig, nrow(res)))
    
    # Dotplot of top terms
    p <- dotplot(res, showCategory = min(15, nrow(res)),
                 title = sprintf("Not-in-FC Tier A genes (n=%d): GO BP", length(tier_a_not_fc))) +
      theme_pub +
      theme(
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        axis.text.y = element_text(size = 7)
      )
  }
  
  add_label(p, "F")
}

# ═══════════════════════════════════════════════════════════════════════════════
# MAIN: Build and Assemble
# ═══════════════════════════════════════════════════════════════════════════════

message("=== Building Figure 6: Functional Coherence ===\n")

message("Panel A: Schematic...")
panel_a <- make_panel_a()

message("Panel B: Tier A up dotplot...")
panel_b <- make_panel_b()

message("Panel C: Tier A down dotplot...")
panel_c <- make_panel_c()

message("Panel D: GSEA curves...")
panel_d <- make_panel_d()

message("Panel E: Coherence barplot...")
panel_e <- make_panel_e()

message("Panel F: Not-in-FC enrichment...")
panel_f <- make_panel_f()

message("\nAssembling figure...")

# Layout: 3 rows × 2 columns
# Row 1: A + B  (schematic + testis dotplot)
# Row 2: C + D  (ovary dotplot + GSEA curves)
# Row 3: E + F  (barplot + not-in-FC enrichment)

design <- "
AABB
CCDD
EEFF
"

figure6 <- wrap_plots(
  A = panel_a, B = panel_b,
  C = panel_c, D = panel_d,
  E = panel_e, F = panel_f,
  design = design
) +
  plot_annotation(
    title = "Figure 6. Functional coherence of consensus tiers in Bombyx testis\u2013ovary transcriptomes",
    theme = theme(
      plot.title = element_text(size = 13, face = "bold", family = "", hjust = 0.5,
                                margin = margin(b = 10))
    )
  )

# ── Save outputs ──────────────────────────────────────────────────────────────

message("\nSaving outputs...")

# Main figure (12×16 inches)
out_main <- file.path(BENCH_FIG, "figure6_functional_coherence.pdf")
cairo_pdf(out_main, width = 12, height = 16)
print(figure6)
dev.off()
message(sprintf("  Main figure: %s", out_main))

# Individual panels
for (nm in c("a", "b", "c", "d", "e", "f")) {
  panel_obj <- get(paste0("panel_", nm))
  out_panel <- file.path(EXP_FIG, sprintf("figure6_panel_%s.pdf", nm))
  cairo_pdf(out_panel, width = 6, height = 5.3)
  print(panel_obj)
  dev.off()
  message(sprintf("  Panel %s: %s", nm, out_panel))
}

# Save schematic data for reproducibility
panel_a_design <- list(
  boxes = data.frame(
    x = c(0.05, 0.38, 0.72, 0.72, 0.72, 0.72),
    y = c(0.55, 0.55, 0.80, 0.58, 0.36, 0.14),
    label = c("Consensus Pipeline", "Tier Classification", "Tier A up (testis-enriched)",
              "Tier A down (ovary-enriched)", "Tier C (low-confidence)", "OFF-only (decontamination-sensitive)"),
    sublabel = c("4 quantifiers: DESeq2 → RRA + CCT", "Evidence-based tier assignment",
                 "n = 2,693; top: cilium assembly", "n = 1,546; top: translation",
                 "n = 507; housekeeping terms", "n = 384; 0 significant terms"),
    fill = c("#333333", "#666666", "#E74C3C", "#3498DB", "#95A5A6", "#BDC3C7")
  )
)
saveRDS(panel_a_design, file.path(EXP_FIG, "figure6_panel_a_design.rds"))
message(sprintf("  Design data: %s", file.path(EXP_FIG, "figure6_panel_a_design.rds")))

message("\nDone. Figure 6 complete.")
