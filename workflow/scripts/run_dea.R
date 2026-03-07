# workflow/scripts/run_dea.R
# Author: 
# Date: 2026-02-09
# Description: Differential Expression Analysis using DESeq2, edgeR, and limma-voom

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(DESeq2)
library(edgeR)
library(limma)
library(readr)
library(tibble)
library(dplyr)
library(tximport)
library(argparse)

# Parameters from Snakemake
input_files <- snakemake@input[["counts"]]
sample_file <- snakemake@params[["sample_file"]]
quant_method <- snakemake@wildcards[["quantifier"]]
gtf_file <- snakemake@input[["gtf"]]
output_dir <- snakemake@output[["outdir"]]
config <- snakemake@config[["dea"]]

# Load sample information
samples <- read_tsv(sample_file, show_col_types = FALSE) %>%
  as.data.frame()
rownames(samples) <- samples$sample

# Determine design formula
group_col <- "group"
if (!group_col %in% colnames(samples)) {
  stop("Column 'group' not found in sample file.")
}

design_formula <- as.formula(paste0("~ ", group_col))
if (!is.null(config$batch_column) && config$batch_column %in% colnames(samples)) {
  design_formula <- as.formula(paste0("~ ", config$batch_column, " + ", group_col))
  cat("Including batch effect in design: ", config$batch_column, "\n")
}

# Load Counts Data
counts_data <- NULL
txi <- NULL

cat("Loading data for method:", quant_method, "\n")

if (quant_method == "featurecounts") {
  # FeatureCounts aggregation usually produces a matrix
  # Format: Geneid, Chr, Start, End, Strand, Length, Sample1, Sample2...
  # We need to skip metadata columns
  # Assuming snakemake input is the matrix file
  
  raw_counts <- read_tsv(input_files, comment = "#", show_col_types = FALSE)
  # Usually featureCounts output has Geneid as first column.
  # The aggregation rule in featurecounts_all outputs "counts_matrix.txt"
  # Let's verify format. Usually featureCounts standard output has headers.
  
  # Standard featureCounts output
  if("Geneid" %in% colnames(raw_counts)){
    counts_data <- raw_counts %>%
      select(-c(Chr, Start, End, Strand, Length)) %>%
      column_to_rownames("Geneid") %>%
      as.matrix()
    
    # Clean up column names for featureCounts (paths to sample names)
    # E.g., results/03.alignment/SAMPLE.bam -> SAMPLE
    colnames(counts_data) <- basename(colnames(counts_data))
    colnames(counts_data) <- gsub("\\.bam$", "", colnames(counts_data))
    
  } else {
     # Maybe the aggregation script simplified it
     counts_data <- raw_counts %>% column_to_rownames(colnames(raw_counts)[1]) %>% as.matrix()
  }

} else if (quant_method == "stringtie") {
  # StringTie aggregation output (gene_counts)
  # Assuming format: gene_id, sample1, sample2...
  raw_counts <- read_tsv(input_files, show_col_types = FALSE)
  counts_data <- raw_counts %>%
    column_to_rownames(colnames(raw_counts)[1]) %>%
    as.matrix()
  
  # Clean up column names for StringTie
  # E.g., SAMPLE/final -> SAMPLE
  colnames(counts_data) <- gsub("/final$", "", colnames(counts_data))

} else if (quant_method %in% c("salmon", "kallisto")) {
  # Use tximport
  files <- input_files # List of files
  names(files) <- samples$sample 
  # Check if files order matches samples
  # The input should be a list of files corresponding to samples.
  # We need a tx2gene mapping.
  
  library(GenomicFeatures)
  # This might be slow. Optimization: Cache or use existing mapping.
  # For now, generate from GTF if not provided.
  # Actually, let's look at the tx2gene file provided in input if any.
  
  if ("tx2gene" %in% names(snakemake@input)) {
     tx2gene <- read_tsv(snakemake@input[["tx2gene"]], col_names=c("tx", "gene"), show_col_types=FALSE)
  } else {
     # Build from GTF
     cat("Building tx2gene from GTF...\n")
     txdb <- makeTxDbFromGFF(gtf_file, format="gtf")
     k <- keys(txdb, keytype = "TXNAME")
     tx2gene <- select(txdb, k, "GENEID", "TXNAME")
  }
  
  type <- ifelse(quant_method == "salmon", "salmon", "kallisto")
  txi <- tximport(files, type = type, tx2gene = tx2gene)
  counts_data <- txi$counts
}

# Ensure columns match samples
common_samples <- intersect(rownames(samples), colnames(counts_data))
if (length(common_samples) == 0) stop("No matching samples between data and config.")
samples <- samples[common_samples, ]
counts_data <- counts_data[, common_samples]

# Keep integers
counts_matrix <- round(counts_data)
mode(counts_matrix) <- "integer"

# Filtering
# Keep genes with at least N counts in at least M samples
# Rule: at least in N samples count >= 10 (N = min group size)
min_group_size <- min(table(samples[[group_col]]))
keep <- rowSums(counts_matrix >= config$min_count) >= min_group_size
counts_matrix <- counts_matrix[keep, ]
cat("Genes after filtering:", nrow(counts_matrix), "\n")

# Run DEA Methods

# Helper to format results
format_res <- function(res, method) {
  res %>%
    as.data.frame() %>%
    rownames_to_column("gene_id") %>%
    mutate(method = method)
}

results_list <- list()

# 1. DESeq2
if ("deseq2" %in% config$methods) {
  cat("Running DESeq2...\n")
  dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                                colData = samples,
                                design = design_formula)
  dds <- DESeq(dds)
  
  # For pairwise comparisons or specific contrasts
  # Simple case: if 2 groups, just results(dds)
  # If >2 groups, integrated result requires contrast handling.
  # User said "pairwise comparison or specified control".
  # For now, let's extract all pairwise comparisons or one default "all".
  # To simplify, we'll iterate pairs if comparisons="all".
  
  groups <- unique(samples[[group_col]])
  
  if (length(groups) >= 2) {
      if (config$comparisons == "all") {
         # Generate all pairs
         pairs <- combn(groups, 2, simplify=FALSE)
      } else {
         # TODO: parse specific comparisons
         # For now assume all pairs
         pairs <- combn(groups, 2, simplify=FALSE)
      }
      
      for (p in pairs) {
          contrast_name <- paste0(p[1], "_vs_", p[2])
          res <- results(dds, contrast = c(group_col, p[1], p[2]))
          res_df <- format_res(res, "deseq2") %>%
              rename(logFC = log2FoldChange, P.Value = pvalue, adj.P.Val = padj) %>%
              mutate(contrast = contrast_name)
          results_list[[paste("deseq2", contrast_name, sep=".")]] <- res_df
      }
  }
}

# 2. edgeR
if ("edger" %in% config$methods) {
  cat("Running edgeR...\n")
  y <- DGEList(counts = counts_matrix, group = samples[[group_col]])
  y <- calcNormFactors(y)
  design <- model.matrix(design_formula, data = samples)
  y <- estimateDisp(y, design)
  
  # Use QL F-test (more robust) or Likelihood Ratio
  fit <- glmQLFit(y, design)
  
  # Contrats
  # Need to construct contrasts manually based on design matrix
  # If simple pairwise, makeContrasts works better with 0-intercept model or just using coefs
  # Lets assume standard design (intercept + factors)
  
  # To handle arbitrary pairs, it's easier to use a 0-intercept model for edgeR/limma
  design_0 <- model.matrix(as.formula(paste0("~ 0 + ", group_col)), data = samples)
  colnames(design_0) <- gsub(group_col, "", colnames(design_0)) # Clean names
  
  # Re-fit with 0-intercept
  y_0 <- estimateDisp(y, design_0) # Re-estimate? Usually minimal specific dispersion change but good practice
  fit_0 <- glmQLFit(y_0, design_0)
  
  if (length(groups) >= 2) {
       for (p in pairs) {
          contrast_vec <- makeContrasts(contrasts = paste0(p[1], "-", p[2]), levels = design_0)
          res <- glmQLFTest(fit_0, contrast = contrast_vec)
          res_tab <- topTags(res, n = Inf)$table
          
          res_df <- res_tab %>%
              rownames_to_column("gene_id") %>%
              mutate(method = "edger") %>%
              rename(logFC = logFC, P.Value = PValue, adj.P.Val = FDR) %>%
              mutate(contrast = paste0(p[1], "_vs_", p[2]))
          
          results_list[[paste("edger", paste0(p[1], "_vs_", p[2]), sep=".")]] <- res_df
       }
  }
}

# 3. limma-voom
if ("limma" %in% config$methods) {
  cat("Running limma-voom...\n")
  dge <- DGEList(counts = counts_matrix)
  dge <- calcNormFactors(dge)
  
  design_0 <- model.matrix(as.formula(paste0("~ 0 + ", group_col)), data = samples)
  colnames(design_0) <- gsub(group_col, "", colnames(design_0))

  v <- voom(dge, design_0, plot = FALSE)
  fit <- lmFit(v, design_0)
  
  if (length(groups) >= 2) {
       for (p in pairs) {
          contrast_vec <- makeContrasts(contrasts = paste0(p[1], "-", p[2]), levels = design_0)
          fit2 <- contrasts.fit(fit, contrast_vec)
          fit2 <- eBayes(fit2)
          res_tab <- topTable(fit2, number = Inf)
          
          res_df <- res_tab %>%
              rownames_to_column("gene_id") %>%
              mutate(method = "limma") %>%
              rename(logFC = logFC, P.Value = P.Value, adj.P.Val = adj.P.Val) %>%
              mutate(contrast = paste0(p[1], "_vs_", p[2]))
          
          results_list[[paste("limma", paste0(p[1], "_vs_", p[2]), sep=".")]] <- res_df
       }
  }
}

# Combine and save
cat("Saving results...\n")
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

for (name in names(results_list)) {
   write_csv(results_list[[name]], file.path(output_dir, paste0(name, ".csv")))
}

# Save normalized counts for visualization
# Using VST from DESeq2 is popular, or logCPM from edgeR
# Let's save VST if DESeq2 was run, else logCPM
if ("deseq2" %in% config$methods) {
  vsd <- vst(dds, blind = FALSE)
  norm_counts <- assay(vsd)
} else {
  norm_counts <- cpm(counts_matrix, log=TRUE)
}

write.csv(norm_counts, file.path(output_dir, "normalized_counts.csv"))
saveRDS(list(counts=counts_matrix, samples=samples, results=results_list, norm_counts=norm_counts), file.path(output_dir, "dea_session.rds"))

sink()
