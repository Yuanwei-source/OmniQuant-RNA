# workflow/scripts/integrate_results.R
# Integrate and Visualize DEA results

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(dplyr)
library(readr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(VennDiagram)

input_rds <- snakemake@input[["rds"]]
output_dir <- snakemake@output[["outdir"]]
config <- snakemake@config[["dea"]]

if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

data <- readRDS(input_rds)
results_list <- data$results
norm_counts <- data$norm_counts
samples <- data$samples

# Extract unique contrasts
contrasts <- unique(sapply(strsplit(names(results_list), "\\."), function(x) paste(x[-1], collapse=".")))

# Thresholds
fdr_th <- config$fdr_threshold
lfc_th <- config$lfc_threshold

for (contrast in contrasts) {
  cat("Processing contrast:", contrast, "\n")
  
  methods <- unique(sapply(strsplit(names(results_list)[grep(contrast, names(results_list))], "\\."), function(x) x[1]))
  
  de_genes_list <- list()
  combined_df <- NULL
  
  for (method in methods) {
    res_name <- paste(method, contrast, sep=".")
    df <- results_list[[res_name]]

    if ("included_in_main" %in% colnames(df)) {
      df <- df %>% filter(included_in_main)
    }
    
    # Filter DE
    # Note: different column names handled in perform_quantifier_dea.R (logFC, P.Value, adj.P.Val)
    # Handle NAs in P-values or logFC by treating them as not DE
    df <- df %>% 
      mutate(is_de = coalesce(adj.P.Val < fdr_th & abs(logFC) > lfc_th, FALSE))
    
    de_genes <- df$gene_id[df$is_de]
    de_genes_list[[method]] <- de_genes
    
    # Save individual DE list
    write_csv(df %>% filter(is_de), file.path(output_dir, paste0("DE_genes_", contrast, "_", method, ".csv")))
    
    # Volcano Plot
    p <- ggplot(df, aes(x=logFC, y=-log10(adj.P.Val), color=is_de)) +
      geom_point(alpha=0.6, size=1) +
      scale_color_manual(values=c("grey", "red")) +
      geom_vline(xintercept=c(-lfc_th, lfc_th), linetype="dashed") +
      geom_hline(yintercept=-log10(fdr_th), linetype="dashed") +
      theme_minimal() +
      labs(title=paste("Volcano Plot:", contrast, method))
    
    # Label top genes
    top_genes <- df %>% filter(is_de) %>% top_n(10, wt=-adj.P.Val)
    p <- p + geom_text_repel(data=top_genes, aes(label=gene_id), max.overlaps=10)
    
    ggsave(file.path(output_dir, paste0("Volcano_", contrast, "_", method, ".pdf")), p)
    
    # Add to combined DF
    sub_df <- df %>% select(gene_id, logFC, adj.P.Val) %>% 
      rename_with(~paste0(., "_", method), -gene_id)
    
    if (is.null(combined_df)) {
      combined_df <- sub_df
    } else {
      combined_df <- full_join(combined_df, sub_df, by="gene_id")
    }
  }
  
  # Venn Diagram
  if (length(methods) > 1 && length(methods) <= 5) {
     pdf(file.path(output_dir, paste0("Venn_", contrast, ".pdf")))
     venn.plot <- venn.diagram(
        x = de_genes_list,
        filename = NULL,
        fill = brewer.pal(length(methods), "Pastel1"),
        disable.logging = TRUE
     )
     grid.draw(venn.plot)
     dev.off()
  }
  
  # Intersection (High Confidence)
  # Genes present in at least min_methods
  all_de_genes <- unlist(de_genes_list)
  gene_counts <- table(all_de_genes)
  high_conf_genes <- names(gene_counts)[gene_counts >= config$min_methods]
  
  write_lines(high_conf_genes, file.path(output_dir, paste0("HighConfidence_DEGs_", contrast, ".txt")))
  
  # Heatmap of High Confidence Genes
  if (length(high_conf_genes) > 1) {
    # Subset normalized counts
    mat <- norm_counts[high_conf_genes, , drop=FALSE]
    
    # Remove rows with zero variance to avoid scale() errors (NaN)
    row_vars <- apply(mat, 1, var)
    keep_rows <- row_vars > 1e-6 & !is.na(row_vars)
    mat <- mat[keep_rows, , drop=FALSE]
    
    if (nrow(mat) > 1) {
        # Scale row
        # Use na.rm=TRUE for scale in case of remaining NAs
        mat_scaled <- t(scale(t(mat)))
        
        # Check for non-finite values (NA, NaN, Inf) and replace
        if(any(!is.finite(mat_scaled))) {
            cat("Warning: Replacing non-finite values in heatmap matrix with 0.\n")
            mat_scaled[!is.finite(mat_scaled)] <- 0
        }
        
        # Explicit check for NAs just in case is.finite missed something (unlikely but safe)
        if(any(is.na(mat_scaled))) {
            cat("Warning: Found NAs in scaled matrix, replacing with 0.\n")
            mat_scaled[is.na(mat_scaled)] <- 0
        }

        # Ensure matrix is numeric
        mode(mat_scaled) <- "numeric"

        cat("Plotting heatmap for", nrow(mat_scaled), "genes.\n")
    
        tryCatch({
            pdf(file.path(output_dir, paste0("Heatmap_", contrast, ".pdf")), width=8, height=10)
            pheatmap(mat_scaled, 
                     annotation_col = samples %>% select(group),
                     show_rownames = length(high_conf_genes) < 50,
                     main = paste("High Confidence DEGs:", contrast),
                     cluster_rows = TRUE,
                     cluster_cols = TRUE)
            dev.off()
        }, error = function(e) {
            cat("Error generating heatmap:", e$message, "\n")
            if(file.exists(file.path(output_dir, paste0("Heatmap_", contrast, ".pdf")))) {
                try(dev.off(), silent=TRUE) 
            }
        })
    }
  }
  
  # Save Combined Results Table
  # Handling case where combined_df might be NULL if no methods ran
  if(!is.null(combined_df)) {
      combined_df <- combined_df %>%
        mutate(HighConfidence = gene_id %in% high_conf_genes)
      
      write_csv(combined_df, file.path(output_dir, paste0("Combined_Results_", contrast, ".csv")))
  }
}

# PCA Plot (Global)
# Clean data for PCA and Correlation
cat("Generating PCA plot...\n")
norm_counts_clean <- norm_counts
norm_counts_clean[!is.finite(norm_counts_clean)] <- 0
# Remove zero variance rows
row_vars_all <- apply(norm_counts_clean, 1, var)
# Handle NAs in row_vars_all which can happen if a row is all NAs
row_vars_all[is.na(row_vars_all)] <- 0
norm_counts_clean <- norm_counts_clean[row_vars_all > 1e-6, ]

if (nrow(norm_counts_clean) > 2) {
    # PCA requires at least as many genes as components usually, but practically >2 is good
    tryCatch({
        pca <- prcomp(t(norm_counts_clean))
        pca_data <- as.data.frame(pca$x)
        pca_data$group <- samples$group
        pca_data$sample <- rownames(pca_data)
        
        p_pca <- ggplot(pca_data, aes(x=PC1, y=PC2, color=group, label=sample)) +
          geom_point(size=3) +
          geom_text_repel() +
          theme_minimal() +
          labs(title="PCA Analysis")
        
        ggsave(file.path(output_dir, "PCA_plot.pdf"), p_pca)
    }, error = function(e) {
         cat("Error generating PCA:", e$message, "\n")
         # Create empty file to avoid snakemake failure
         file.create(file.path(output_dir, "PCA_plot.pdf"))
    })
    
    # Sample Correlation
    tryCatch({
        sample_dists <- dist(t(norm_counts_clean))
        sample_dist_matrix <- as.matrix(sample_dists)
        
        pdf(file.path(output_dir, "Sample_Correlation.pdf"))
        # Check if we have enough variability/valid values
        if(all(!is.na(sample_dist_matrix)) && all(is.finite(sample_dist_matrix))) {
            pheatmap(sample_dist_matrix, clustering_distance_rows=sample_dists, clustering_distance_cols=sample_dists, main="Sample-to-Sample Distances")
        } else {
             cat("Skipping sample correlation heatmap due to NA/Inf values in distance matrix.\n")
             plot(1, type="n", axes=FALSE, xlab="", ylab="", main="Correlation Failed: NAs in distance matrix")
        }
        dev.off()
    }, error = function(e) {
         cat("Error generating Correlation plot:", e$message, "\n")
         if(file.exists(file.path(output_dir, "Sample_Correlation.pdf"))) { try(dev.off(), silent=TRUE) }
    })

} else {
    cat("No variable genes found for PCA/Correlation.\n")
    # Touch output to avoid Snakemake error
    file.create(file.path(output_dir, "PCA_plot.pdf"))
}

sink()
