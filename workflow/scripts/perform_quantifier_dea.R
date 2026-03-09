# workflow/scripts/perform_quantifier_dea.R
# Unified differential expression analysis for a single quantifier using featureCounts or tximport inputs

log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)
  library(limma)
  library(readr)
  library(tibble)
  library(dplyr)
  library(tximport)
})

invisible(utils::globalVariables(c(
  "gene_id", "gene_name", "is_reference_gene", "allow_consensus_main",
  "is_novel_only", "is_ambiguous", "gene_id_original", "gene_id_resolved",
  "transcript_id", "gene_id_standard", "resolution_status",
  "gene_name_namespace", "is_reference_gene_namespace",
  "included_in_main_namespace", "included_in_main"
)))

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

make_import_summary <- function(metric, value) {
  tibble(
    metric = as.character(metric),
    value = as.character(value)
  )
}

as_flag <- function(x) {
  if (is.logical(x)) return(isTRUE(x))
  if (is.null(x) || length(x) == 0) return(FALSE)
  tolower(as.character(x)[1]) %in% c("true", "1", "yes")
}

load_samples <- function(sample_file) {
  samples <- read_tsv(sample_file, show_col_types = FALSE) %>% as.data.frame()
  if ("group" %in% colnames(samples)) {
    samples$group <- as.factor(samples$group)
  }
  rownames(samples) <- samples$sample
  samples
}

build_pairs <- function(groups, comparisons) {
  groups <- as.character(groups)
  if (length(groups) < 2) {
    stop("At least two groups are required for DEA.")
  }

  if (is.character(comparisons) && length(comparisons) == 1 && comparisons == "all") {
    return(combn(groups, 2, simplify = FALSE))
  }

  pairs <- list()
  comparison_values <- comparisons
  if (!is.vector(comparison_values)) {
    comparison_values <- as.vector(comparison_values)
  }

  for (comp in comparison_values) {
    parts <- strsplit(as.character(comp), "_vs_")[[1]]
    if (length(parts) != 2) {
      warning("Skipping malformed comparison: ", comp)
      next
    }
    if (!all(parts %in% groups)) {
      warning("Skipping comparison with unknown groups: ", comp)
      next
    }
    pairs[[length(pairs) + 1]] <- parts
  }

  if (length(pairs) == 0) {
    stop("No valid comparisons could be constructed from configuration.")
  }
  pairs
}

load_gene_namespace <- function(namespace_path) {
  namespace <- read_tsv(namespace_path, show_col_types = FALSE) %>%
    transmute(
      gene_id_standard = .data$gene_id,
      gene_name = coalesce(na_if(.data$gene_name, "."), .data$gene_id),
      is_reference_gene = as.logical(toupper(.data$is_reference_gene) == "TRUE"),
      included_in_main = as.logical(toupper(.data$allow_consensus_main) == "TRUE")
    )
  namespace
}

load_featurecounts_matrix <- function(counts_path, gene_namespace) {
  raw_counts <- read_tsv(counts_path, comment = "#", show_col_types = FALSE)

  if ("Geneid" %in% colnames(raw_counts)) {
    counts_data <- raw_counts %>%
      select(-any_of(c("Chr", "Start", "End", "Strand", "Length"))) %>%
      column_to_rownames("Geneid") %>%
      as.matrix()
  } else {
    counts_data <- raw_counts %>% column_to_rownames(colnames(raw_counts)[1]) %>% as.matrix()
  }

  colnames(counts_data) <- basename(colnames(counts_data))
  colnames(counts_data) <- gsub("\\.bam$", "", colnames(counts_data))

  gene_metadata <- gene_namespace %>%
    mutate(
      is_novel_only = FALSE,
      is_ambiguous = FALSE,
      mapping_policy = "raw_gene_counts"
    )

  import_summary <- make_import_summary(
    metric = c("input_mode", "features_loaded", "reference_genes_loaded"),
    value = c("gene_counts_matrix", nrow(counts_data), sum(gene_metadata$is_reference_gene))
  )

  list(
    counts = counts_data,
    txi_raw = NULL,
    txi_scaled = NULL,
    gene_metadata = gene_metadata,
    import_summary = import_summary
  )
}

load_tx2gene_master <- function(master_path, quantifier, main_only = TRUE) {
  master <- read_tsv(master_path, show_col_types = FALSE) %>%
    mutate(
      is_reference_gene = as.logical(toupper(.data$is_reference_gene) == "TRUE"),
      is_novel_only = as.logical(toupper(.data$is_novel_only) == "TRUE"),
      is_ambiguous = as.logical(toupper(.data$is_ambiguous) == "TRUE"),
      allow_consensus_main = as.logical(toupper(.data$allow_consensus_main) == "TRUE"),
      gene_name = coalesce(na_if(.data$gene_name, "."), .data$gene_id_original, .data$gene_id_resolved)
    )

  quantifier_key <- if (quantifier %in% c("salmon", "kallisto")) "reference" else "stringtie"
  quant_master <- master %>% filter(.data$quantifier == quantifier_key)
  filtered_master <- quant_master
  if (main_only) {
    filtered_master <- filtered_master %>% filter(.data$allow_consensus_main)
  }

  tx2gene <- filtered_master %>%
    filter(!is.na(.data$gene_id_resolved), .data$gene_id_resolved != "") %>%
    distinct(.data$transcript_id, .data$gene_id_resolved)

  gene_metadata <- quant_master %>%
    mutate(gene_id_standard = coalesce(.data$gene_id_resolved, .data$gene_id_original)) %>%
    distinct(.data$gene_id_standard, .keep_all = TRUE) %>%
    transmute(
      gene_id_standard = .data$gene_id_standard,
      gene_name = coalesce(na_if(.data$gene_name, "."), .data$gene_id_standard),
      is_reference_gene = .data$is_reference_gene,
      is_novel_only = .data$is_novel_only,
      is_ambiguous = .data$is_ambiguous,
      included_in_main = .data$allow_consensus_main,
      mapping_policy = .data$resolution_status
    )

  import_summary <- make_import_summary(
    metric = c(
      "master_rows_total",
      "master_rows_filtered",
      "tx2gene_rows_used",
      "novel_only_rows",
      "ambiguous_rows"
    ),
    value = c(
      nrow(quant_master),
      nrow(filtered_master),
      nrow(tx2gene),
      sum(quant_master$is_novel_only),
      sum(quant_master$is_ambiguous)
    )
  )

  list(
    tx2gene = tx2gene,
    gene_metadata = gene_metadata,
    import_summary = import_summary
  )
}

internal_import_function <- function(input_mode, primary_input, tx2gene_master_path, gene_namespace_path,
                                     quantifier, read_length, main_only) {
  gene_namespace <- load_gene_namespace(gene_namespace_path)

  if (input_mode == "gene_counts_matrix") {
    return(load_featurecounts_matrix(primary_input, gene_namespace))
  }

  manifest <- read_tsv(primary_input, show_col_types = FALSE)
  files <- manifest$file_path
  names(files) <- manifest$sample

  master_info <- load_tx2gene_master(tx2gene_master_path, quantifier, main_only = main_only)
  tx2gene <- as.data.frame(master_info$tx2gene)
  colnames(tx2gene) <- c("TXNAME", "GENEID")

  import_type <- unique(manifest$import_type)
  if (length(import_type) != 1) {
    stop("Manifest must contain exactly one import_type. Got: ", paste(import_type, collapse = ", "))
  }

  tximport_args <- list(files = files, type = import_type, tx2gene = tx2gene)
  if (identical(import_type, "stringtie")) {
    tximport_args$readLength <- read_length
  }
  if (identical(import_type, "kallisto")) {
    tximport_args$ignoreAfterBar <- TRUE
  }

  txi_raw <- do.call(tximport::tximport, c(tximport_args, list(countsFromAbundance = "no")))
  txi_scaled <- do.call(tximport::tximport, c(tximport_args, list(countsFromAbundance = "lengthScaledTPM")))

  gene_metadata <- master_info$gene_metadata %>%
    full_join(gene_namespace, by = "gene_id_standard", suffix = c("", "_namespace")) %>%
    transmute(
      gene_id_standard = .data$gene_id_standard,
      gene_name = coalesce(.data$gene_name_namespace, .data$gene_name, .data$gene_id_standard),
      is_reference_gene = coalesce(.data$is_reference_gene_namespace, .data$is_reference_gene, FALSE),
      is_novel_only = coalesce(.data$is_novel_only, FALSE),
      is_ambiguous = coalesce(.data$is_ambiguous, FALSE),
      included_in_main = coalesce(.data$included_in_main, .data$included_in_main_namespace, FALSE),
      mapping_policy = coalesce(.data$mapping_policy, "tximport")
    )

  import_summary <- bind_rows(
    make_import_summary(
      metric = c("input_mode", "samples_loaded", "genes_imported"),
      value = c(input_mode, length(files), nrow(txi_raw$counts))
    ),
    master_info$import_summary
  )

  list(
    counts = txi_raw$counts,
    txi_raw = txi_raw,
    txi_scaled = txi_scaled,
    gene_metadata = gene_metadata,
    import_summary = import_summary
  )
}

build_tximport_edgeR <- function(txi) {
  cts <- txi$counts
  normMat <- txi$length
  normMat <- normMat / exp(rowMeans(log(normMat)))
  normCts <- cts / normMat
  eff.lib <- calcNormFactors(normCts) * colSums(normCts)
  normMat <- sweep(normMat, 2, eff.lib, "*")
  normMat <- log(normMat)
  y <- DGEList(cts)
  y <- scaleOffset(y, normMat)
  y
}

annotate_result <- function(df, gene_metadata, quantifier, method, contrast, input_mode, mapping_policy, main_only) {
  metadata <- gene_metadata %>% distinct(.data$gene_id_standard, .keep_all = TRUE)
  df %>%
    mutate(gene_id_standard = .data$gene_id) %>%
    left_join(metadata, by = "gene_id_standard") %>%
    mutate(
      gene_name = coalesce(.data$gene_name, .data$gene_id_standard),
      quantifier = quantifier,
      method = method,
      contrast = contrast,
      input_mode = input_mode,
      is_reference_gene = coalesce(.data$is_reference_gene, FALSE),
      is_novel_only = coalesce(.data$is_novel_only, FALSE),
      is_ambiguous = coalesce(.data$is_ambiguous, FALSE),
      included_in_main = if (main_only) TRUE else coalesce(.data$included_in_main, FALSE),
      mapping_policy = mapping_policy
    ) %>%
    select(
      all_of(c("gene_id", "gene_id_standard", "gene_name")),
      everything()
    )
}

# Parameters from Snakemake
primary_input <- snakemake@input[["primary"]]
sample_file <- snakemake@input[["sample_file"]]
tx2gene_master_path <- snakemake@input[["tx2gene_master"]]
gene_namespace_path <- snakemake@input[["gene_namespace"]]
quant_method <- snakemake@wildcards[["quantifier"]]
output_dir <- snakemake@output[["outdir"]]
normalized_counts_output <- snakemake@output[["norm_counts"]]
import_summary_output <- snakemake@output[["import_summary"]]
config_dea <- snakemake@config[["dea"]]
input_mode <- snakemake@params[["input_mode"]]
read_length <- as.numeric(snakemake@params[["read_length"]])
main_only <- as_flag(snakemake@params[["main_only"]])
mapping_policy <- snakemake@params[["mapping_policy"]]
novel_policy <- snakemake@params[["novel_policy"]]

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

samples <- load_samples(sample_file)
group_col <- "group"
if (!group_col %in% colnames(samples)) {
  stop("Column 'group' not found in sample file.")
}

design_formula <- as.formula(paste0("~ ", group_col))
if (!is.null(config_dea$batch_column) && config_dea$batch_column %in% colnames(samples)) {
  design_formula <- as.formula(paste0("~ ", config_dea$batch_column, " + ", group_col))
  cat("Including batch effect in design: ", config_dea$batch_column, "\n")
}

cat("Loading data for method:", quant_method, "\n")
imported <- internal_import_function(
  input_mode = input_mode,
  primary_input = primary_input,
  tx2gene_master_path = tx2gene_master_path,
  gene_namespace_path = gene_namespace_path,
  quantifier = quant_method,
  read_length = read_length,
  main_only = main_only
)

counts_data <- imported$counts
common_samples <- intersect(rownames(samples), colnames(counts_data))
if (length(common_samples) == 0) stop("No matching samples between data and config.")
samples <- samples[common_samples, , drop = FALSE]
counts_data <- counts_data[, common_samples, drop = FALSE]
if (!is.null(imported$txi_raw)) {
  imported$txi_raw$counts <- imported$txi_raw$counts[, common_samples, drop = FALSE]
  imported$txi_raw$abundance <- imported$txi_raw$abundance[, common_samples, drop = FALSE]
  imported$txi_raw$length <- imported$txi_raw$length[, common_samples, drop = FALSE]
  imported$txi_scaled$counts <- imported$txi_scaled$counts[, common_samples, drop = FALSE]
  imported$txi_scaled$abundance <- imported$txi_scaled$abundance[, common_samples, drop = FALSE]
  imported$txi_scaled$length <- imported$txi_scaled$length[, common_samples, drop = FALSE]
}

counts_matrix <- round(counts_data)
mode(counts_matrix) <- "integer"

gene_metadata <- imported$gene_metadata %>% distinct(gene_id_standard, .keep_all = TRUE)

min_group_size <- min(table(samples[[group_col]]))
keep <- rowSums(counts_matrix >= config_dea$min_count) >= min_group_size
counts_matrix <- counts_matrix[keep, , drop = FALSE]
gene_metadata <- gene_metadata %>% filter(gene_id_standard %in% rownames(counts_matrix))
cat("Genes after filtering:", nrow(counts_matrix), "\n")

write_tsv(
  bind_rows(imported$import_summary, make_import_summary(metric = "genes_after_filtering", value = nrow(counts_matrix))),
  import_summary_output
)

results_list <- list()
groups <- as.character(unique(samples[[group_col]]))
pairs <- build_pairs(groups, config_dea$comparisons)

format_res <- function(res, method, contrast_name) {
  annotate_result(
    res %>% as.data.frame() %>% rownames_to_column("gene_id"),
    gene_metadata = gene_metadata,
    quantifier = quant_method,
    method = method,
    contrast = contrast_name,
    input_mode = input_mode,
    mapping_policy = mapping_policy,
    main_only = main_only
  )
}

# 1. DESeq2
if ("deseq2" %in% config_dea$methods) {
  cat("Running DESeq2...\n")
  if (!is.null(imported$txi_raw)) {
    txi_raw <- imported$txi_raw
    txi_raw$counts <- txi_raw$counts[rownames(counts_matrix), common_samples, drop = FALSE]
    txi_raw$abundance <- txi_raw$abundance[rownames(counts_matrix), common_samples, drop = FALSE]
    txi_raw$length <- txi_raw$length[rownames(counts_matrix), common_samples, drop = FALSE]
    dds <- DESeqDataSetFromTximport(txi = txi_raw, colData = samples, design = design_formula)
  } else {
    dds <- DESeqDataSetFromMatrix(countData = counts_matrix, colData = samples, design = design_formula)
  }

  dds <- DESeq(dds)
  for (p in pairs) {
    contrast_name <- paste0(p[1], "_vs_", p[2])
    res <- results(dds, contrast = c(group_col, p[1], p[2]))
    res_df <- format_res(res, "deseq2", contrast_name) %>%
      rename(logFC = log2FoldChange, P.Value = pvalue, adj.P.Val = padj)
    results_list[[paste("deseq2", contrast_name, sep = ".")]] <- res_df
  }
}

# 2. edgeR
if ("edger" %in% config_dea$methods) {
  cat("Running edgeR...\n")
  design_0 <- model.matrix(as.formula(paste0("~ 0 + ", group_col)), data = samples)
  colnames(design_0) <- gsub(group_col, "", colnames(design_0))

  if (!is.null(imported$txi_raw)) {
    y <- build_tximport_edgeR(imported$txi_raw)
    y <- y[rownames(counts_matrix), , keep.lib.sizes = FALSE]
  } else {
    y <- DGEList(counts = counts_matrix, group = samples[[group_col]])
    y <- calcNormFactors(y)
  }

  y <- estimateDisp(y, design_0)
  fit_0 <- glmQLFit(y, design_0)

  for (p in pairs) {
    contrast_name <- paste0(p[1], "_vs_", p[2])
    contrast_vec <- makeContrasts(contrasts = paste0(p[1], "-", p[2]), levels = design_0)
    res <- glmQLFTest(fit_0, contrast = contrast_vec)
    res_tab <- topTags(res, n = Inf)$table
    res_df <- annotate_result(
      res_tab %>% rownames_to_column("gene_id") %>% rename(P.Value = PValue, adj.P.Val = FDR),
      gene_metadata = gene_metadata,
      quantifier = quant_method,
      method = "edger",
      contrast = contrast_name,
      input_mode = input_mode,
      mapping_policy = mapping_policy,
      main_only = main_only
    )
    results_list[[paste("edger", contrast_name, sep = ".")]] <- res_df
  }
}

# 3. limma-voom
if ("limma" %in% config_dea$methods) {
  cat("Running limma-voom...\n")
  design_0 <- model.matrix(as.formula(paste0("~ 0 + ", group_col)), data = samples)
  colnames(design_0) <- gsub(group_col, "", colnames(design_0))

  limma_counts <- if (!is.null(imported$txi_scaled)) round(imported$txi_scaled$counts[rownames(counts_matrix), common_samples, drop = FALSE]) else counts_matrix
  dge <- DGEList(counts = limma_counts)
  keep_limma <- rownames(dge$counts) %in% rownames(counts_matrix)
  dge <- dge[keep_limma, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  v <- voom(dge, design_0, plot = FALSE)
  fit <- lmFit(v, design_0)

  for (p in pairs) {
    contrast_name <- paste0(p[1], "_vs_", p[2])
    contrast_vec <- makeContrasts(contrasts = paste0(p[1], "-", p[2]), levels = design_0)
    fit2 <- contrasts.fit(fit, contrast_vec)
    fit2 <- eBayes(fit2)
    res_tab <- topTable(fit2, number = Inf)
    res_df <- annotate_result(
      res_tab %>% rownames_to_column("gene_id"),
      gene_metadata = gene_metadata,
      quantifier = quant_method,
      method = "limma",
      contrast = contrast_name,
      input_mode = input_mode,
      mapping_policy = mapping_policy,
      main_only = main_only
    )
    results_list[[paste("limma", contrast_name, sep = ".")]] <- res_df
  }
}

cat("Saving results...\n")
for (name in names(results_list)) {
  write_csv(results_list[[name]], file.path(output_dir, paste0(name, ".csv")))
}

if ("deseq2" %in% config_dea$methods) {
  vsd <- vst(dds, blind = FALSE)
  norm_counts <- assay(vsd)
} else if (!is.null(imported$txi_scaled)) {
  norm_counts <- imported$txi_scaled$counts[rownames(counts_matrix), common_samples, drop = FALSE]
} else {
  norm_counts <- cpm(counts_matrix, log = TRUE)
}

write.csv(norm_counts, normalized_counts_output)
saveRDS(
  list(
    counts = counts_matrix,
    samples = samples,
    results = results_list,
    norm_counts = norm_counts,
    input_mode = input_mode,
    mapping_policy = mapping_policy,
    novel_policy = novel_policy,
    import_summary = read_tsv(import_summary_output, show_col_types = FALSE)
  ),
  file.path(output_dir, "dea_session.rds")
)

sink()
