#!/usr/bin/env Rscript
###############################################################################
# OmniQuant-RNA: Figure 6 Reassembly — Functional Coherence
#
# Reassembles all 6 panels (A-F) into the final figure.
# Panels A, B, C, E are from existing files.
# Panels D and F are from newly generated files (updated combined GO).
#
# Layout: 3 rows × 2 columns
#   Row 1: A | B
#   Row 2: C | D
#   Row 3: E | F
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(pdftools)
  library(cowplot)
})

# ═════════════════════════════════════════════════════════════════════════════
# 0. Configuration
# ═════════════════════════════════════════════════════════════════════════════

BASE      <- Sys.getenv("PROJECT_ROOT", normalizePath("."))
FIG_DIR   <- file.path(BASE, "experiments/bombyx_enrichment/results/figures")
OUT_DIR   <- file.path(BASE, "benchmark_results/figures")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

OUTPUT <- file.path(OUT_DIR, "figure6_functional_coherence.pdf")

# Panel paths
panels <- list(
  A = file.path(FIG_DIR, "figure6_panel_a.pdf"),
  B = file.path(FIG_DIR, "figure6_panel_b.pdf"),
  C = file.path(FIG_DIR, "figure6_panel_c.pdf"),
  D = file.path(FIG_DIR, "figure6_panel_d.pdf"),
  E = file.path(FIG_DIR, "figure6_panel_e.pdf"),
  F = file.path(FIG_DIR, "figure6_panel_f.pdf")
)

# Verify all panels exist
for (nm in names(panels)) {
  if (!file.exists(panels[[nm]])) {
    stop(sprintf("Panel %s not found: %s", nm, panels[[nm]]))
  }
}
message("All 6 panels found.")

# ═════════════════════════════════════════════════════════════════════════════
# 1. Render each PDF panel as a raster grob for cowplot
# ═════════════════════════════════════════════════════════════════════════════

message("Rendering panels...")

render_panel <- function(pdf_path, dpi = 200) {
  # Render first page of PDF to bitmap, convert to grob
  img <- pdftools::pdf_render_page(pdf_path, page = 1, dpi = dpi)
  cowplot::ggdraw() + cowplot::draw_image(img)
}

panel_grobs <- lapply(panels, render_panel)

# ═════════════════════════════════════════════════════════════════════════════
# 2. Assemble 3×2 layout
# ═════════════════════════════════════════════════════════════════════════════

message("Assembling 3×2 layout...")

# Row 1: A | B
row1 <- cowplot::plot_grid(
  panel_grobs[["A"]], panel_grobs[["B"]],
  ncol = 2,
  labels = c("A", "B"),
  label_size = 14,
  label_fontface = "bold"
)

# Row 2: C | D
row2 <- cowplot::plot_grid(
  panel_grobs[["C"]], panel_grobs[["D"]],
  ncol = 2,
  labels = c("C", "D"),
  label_size = 14,
  label_fontface = "bold"
)

# Row 3: E | F
row3 <- cowplot::plot_grid(
  panel_grobs[["E"]], panel_grobs[["F"]],
  ncol = 2,
  labels = c("E", "F"),
  label_size = 14,
  label_fontface = "bold"
)

# Stack rows
figure6 <- cowplot::plot_grid(
  row1, row2, row3,
  ncol = 1,
  rel_heights = c(1, 1, 1)
)

# ═════════════════════════════════════════════════════════════════════════════
# 3. Save
# ═════════════════════════════════════════════════════════════════════════════

message(sprintf("Saving to: %s", OUTPUT))
pdf(OUTPUT, width = 12, height = 16)
print(figure6)
dev.off()

# Also save PNG for markdown/inline viewing
png_output <- sub("\\.pdf$", ".png", OUTPUT)
png(png_output, width = 12, height = 16, units = "in", res = 150)
print(figure6)
dev.off()
message(sprintf("Saved PNG: %s", png_output))

message("=== Figure 6 reassembly COMPLETE ===")
message(sprintf("  Output: %s", OUTPUT))
