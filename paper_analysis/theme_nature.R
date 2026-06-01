###############################################################################
# Nature publication theme for OmniQuant-RNA figures
# Source this file in any plotting script to apply the unified Nature style.
#
# Usage:
#   source("theme_nature.R")          # from same directory
#   source("../theme_nature.R")       # from subdirectory
#
# Then use save_pub_r(plot, "output_name") to export SVG/PDF/TIFF together.
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
})

# ── Unified palette ───────────────────────────────────────────────────────────
palette_nature <- c(
  neutral_dark   = "#272727",
  neutral_mid    = "#767676",
  neutral_light  = "#D8D8D8",
  signal_blue    = "#3182BD",
  signal_teal    = "#33B5A5",
  accent_red     = "#D24B40",
  accent_orange  = "#E28E2C",
  accent_purple  = "#8C5AA6",
  accent_green   = "#5DA85C"
)

# Convenience aliases for common plot roles
col_consensus  <- palette_nature[["accent_red"]]
col_featurecounts <- palette_nature[["signal_blue"]]
col_salmon     <- palette_nature[["accent_orange"]]
col_kallisto   <- palette_nature[["accent_purple"]]
col_stringtie  <- palette_nature[["signal_teal"]]

col_testis     <- "#3B7DD8"
col_ovary      <- "#E24B4B"
col_tiera      <- "#E41A1C"
col_tierb      <- "#377EB8"
col_tierc      <- "#4DAF4A"

# ── Nature theme ──────────────────────────────────────────────────────────────
theme_nature <- function(base_size = 6.5, base_family = "") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line        = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks       = element_line(linewidth = 0.35, colour = "black"),
      axis.title       = element_text(size = base_size),
      axis.text        = element_text(size = base_size - 0.5),
      legend.title     = element_text(size = base_size - 0.3),
      legend.text      = element_text(size = base_size - 0.7),
      legend.key.size  = unit(0.4, "cm"),
      strip.text       = element_text(size = base_size - 0.3, face = "bold"),
      plot.title       = element_text(size = base_size + 1, face = "bold"),
      plot.subtitle    = element_text(size = base_size - 1),
      plot.caption     = element_text(size = base_size - 2, hjust = 1),
      panel.grid       = element_blank(),
      panel.spacing    = unit(1.2, "mm"),
      plot.margin      = margin(4, 4, 4, 4)
    )
}

# Set as default for all subsequent ggplot calls
theme_set(theme_nature())

# ── Publication export ────────────────────────────────────────────────────────
save_pub_r <- function(plot, filename, width_mm = 183, height_mm = 120, dpi = 600) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4

  # SVG (editable text, vector)
  if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(paste0(filename, ".svg"), width = w, height = h)
    print(plot)
    dev.off()
  }

  # PDF (cairo, editable TrueType)
  grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, family = "Helvetica")
  print(plot)
  dev.off()

  # High-DPI TIFF (raster, for submission)
  has_ragg <- requireNamespace("ragg", quietly = TRUE)
  if (has_ragg) {
    ragg::agg_tiff(paste0(filename, ".tiff"), width = w, height = h,
                    units = "in", res = dpi)
    print(plot)
    dev.off()
  }

  # PNG for quick preview
  if (has_ragg) {
    ragg::agg_png(paste0(filename, ".png"), width = w, height = h,
                   units = "in", res = 150)
    print(plot)
    dev.off()
  } else {
    grDevices::png(paste0(filename, ".png"), width = w, height = h,
                    units = "in", res = 150)
    print(plot)
    dev.off()
  }

  message(sprintf("  Exported: %s.{svg,pdf,tiff,png}", filename))
}

# ── Panel label helper ───────────────────────────────────────────────────────
add_panel_labels <- function(plot, labels = "a", size = 8) {
  plot + plot_annotation(tag_levels = labels) &
    theme(plot.tag = element_text(size = size, face = "bold"))
}

message("Nature theme loaded. Use theme_nature() + save_pub_r() for figures.")
