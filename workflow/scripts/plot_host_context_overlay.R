log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(readr)
})

burden_path <- snakemake@input[["burden"]]
targets_path <- snakemake@input[["targets"]]
sample_path <- snakemake@input[["sample_file"]]
counts_path <- snakemake@input[["normalized_counts"]]
output_path <- snakemake@output[["plot"]]

burden <- read_tsv(burden_path, show_col_types = FALSE)
targets <- read_tsv(targets_path, show_col_types = FALSE)
samples <- read_tsv(sample_path, show_col_types = FALSE)
counts <- read.csv(counts_path, row.names = 1, check.names = FALSE)

burden$non_host_fraction <- as.numeric(burden$non_host_fraction)
targets$virus_presence <- targets$virus_presence == "yes"

shared_samples <- intersect(colnames(counts), burden$sample)
counts <- counts[, shared_samples, drop = FALSE]
burden <- burden[match(shared_samples, burden$sample), , drop = FALSE]
targets <- targets[match(shared_samples, targets$sample), , drop = FALSE]
samples <- samples[match(shared_samples, samples$sample), , drop = FALSE]

row_vars <- apply(counts, 1, var, na.rm = TRUE)
row_vars[is.na(row_vars)] <- 0
counts <- counts[row_vars > 1e-6, , drop = FALSE]

if (nrow(counts) < 2 || ncol(counts) < 2) {
  p <- ggplot() +
    annotate("text", x = 0, y = 0, label = "Insufficient data for host context PCA") +
    theme_void()
  ggsave(output_path, p, width = 8, height = 5)
} else {
  pca <- prcomp(t(counts), scale. = FALSE)
  pca_data <- data.frame(
    sample = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    stringsAsFactors = FALSE
  )

  pca_data <- merge(pca_data, burden, by = "sample", all.x = TRUE)
  pca_data <- merge(pca_data, targets[, c("sample", "virus_presence")], by = "sample", all.x = TRUE)
  pca_data$group <- as.factor(pca_data$group)
  pca_data$virus_presence[is.na(pca_data$virus_presence)] <- FALSE

  p <- ggplot(
    pca_data,
    aes(
      x = PC1,
      y = PC2,
      color = group,
      size = non_host_fraction,
      shape = virus_presence,
      label = sample
    )
  ) +
    geom_point(alpha = 0.9) +
    geom_text_repel(max.overlaps = 20) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17)) +
    scale_size_continuous(range = c(2.5, 8)) +
    labs(
      title = "Host PCA With Microbial Burden Overlay",
      x = "PC1",
      y = "PC2",
      color = "Group",
      size = "Non-host fraction",
      shape = "Virus signal"
    ) +
    theme_minimal(base_size = 11)

  ggsave(output_path, p, width = 8, height = 6)
}

sink(type = "message")
sink()
close(log)