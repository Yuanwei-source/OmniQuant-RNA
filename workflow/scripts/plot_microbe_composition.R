log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
})

input_path <- snakemake@input[["composition"]]
output_path <- snakemake@output[["plot"]]

df <- read_tsv(input_path, show_col_types = FALSE)
df$sample <- factor(df$sample, levels = unique(df$sample))
df$category <- factor(
  df$category,
  levels = c("Technical", "Bacteria", "Fungi", "Viruses", "Other_Classified", "Unclassified")
)
df$fraction_total <- as.numeric(df$fraction_total)

p <- ggplot(df, aes(x = sample, y = fraction_total, fill = category)) +
  geom_col(width = 0.8) +
  scale_fill_manual(
    values = c(
      Technical = "#8c8c8c",
      Bacteria = "#1b9e77",
      Fungi = "#d95f02",
      Viruses = "#7570b3",
      Other_Classified = "#66a61e",
      Unclassified = "#e6ab02"
    )
  ) +
  labs(
    title = "Microbial Clue Composition Across Samples",
    x = "Sample",
    y = "Fraction of Total Read Pairs",
    fill = "Category"
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(output_path, p, width = 9, height = 5)

sink(type = "message")
sink()
close(log)