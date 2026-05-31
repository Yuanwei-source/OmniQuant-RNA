#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(RobustRankAggreg); library(dplyr); library(readr) })

BASE <- Sys.getenv("PROJECT_ROOT", normalizePath("."))
CONSENSUS <- file.path(BASE,
  "benchmark_results/bombyx_mori/runs/2026-05-24_bombyx_mori_decontam-on",
  "07.consensus_expression/Testis_vs_Ovary/consensus_results.tsv")
N_PERM <- 1000L; ALPHA <- 0.05; P_CLIP <- 1e-16

consensus <- read_tsv(CONSENSUS, col_types=cols(.default=col_character()))
pcols <- c("P.Value__featurecounts","P.Value__stringtie","P.Value__salmon","P.Value__kallisto")
lfcols <- c("logFC__featurecounts","logFC__stringtie","logFC__salmon","logFC__kallisto")
for (col in c(pcols, lfcols)) consensus[[col]] <- as.numeric(consensus[[col]])
valid <- complete.cases(consensus[, c(pcols, lfcols)])
data <- consensus[valid, ]; n_genes <- nrow(data)
gene_ids <- data$gene_id_standard
p_mat <- as.matrix(data[, pcols]); lf_mat <- as.matrix(data[, lfcols])

# Pipeline's exact CCT logic (from run_consensus_dea.R compute_cct_scores)
run_cct <- function(p_matrix, lf_matrix, direction, eps=1e-16) {
  mask <- !is.na(p_matrix) & !is.na(lf_matrix)
  if (direction == "up") mask <- mask & (lf_matrix > 0)
  else mask <- mask & (lf_matrix < 0)
  valid_p <- p_matrix
  valid_p[!mask] <- NA
  valid_p[valid_p < eps] <- eps
  tan_mat <- tan((0.5 - valid_p) * pi)
  t_stat <- rowMeans(tan_mat, na.rm = TRUE)
  t_stat[!is.finite(t_stat)] <- 0
  cct_p <- 0.5 - atan(t_stat) / pi
  cct_p[is.nan(cct_p)] <- 1
  pmax(pmin(cct_p, 1), eps)
}

set.seed(42)
rra_sig <- integer(N_PERM); cct_sig <- integer(N_PERM); dual_sig <- integer(N_PERM)
t0 <- Sys.time()

for (i in seq_len(N_PERM)) {
  # Shuffle P-values independently per quantifier
  p_shuf <- p_mat
  for (q in 1:4) p_shuf[, q] <- p_mat[sample(n_genes), q]
  
  # RRA (using shuffled P-values; keep original logFC for direction)
  # Build rank lists: for each quantifier, rank genes by shuffled P-value
  rra_scores_all <- c()
  for (direction in c("up","down")) {
    rank_lists <- vector("list", 4L)
    for (q in 1:4) {
      # One-tailed P from shuffled data
      p_one <- if(direction=="up") ifelse(lf_mat[,q]>0, p_shuf[,q]/2, 1-p_shuf[,q]/2)
              else ifelse(lf_mat[,q]<0, p_shuf[,q]/2, 1-p_shuf[,q]/2)
      ord <- order(p_one)
      rank_lists[[q]] <- gene_ids[ord]
    }
    rra <- tryCatch(aggregateRanks(rank_lists, method="RRA", N=n_genes), error=function(e)NULL)
    if (!is.null(rra) && nrow(rra)>0) rra_scores_all <- c(rra_scores_all, rra$Score)
  }
  rra_sig[i] <- sum(rra_scores_all < ALPHA)
  
  # CCT: use pipeline's logic directly, with shuffled P-values
  cct_up <- run_cct(p_shuf, lf_mat, "up", P_CLIP)
  cct_down <- run_cct(p_shuf, lf_mat, "down", P_CLIP)
  # best_cct = min of up/down per gene
  cct_p <- pmin(cct_up, cct_down)
  cct_fdr <- p.adjust(cct_p, method="BH")
  cct_sig[i] <- sum(cct_fdr < ALPHA)
  
  # Dual: must pass both RRA AND CCT
  # For each gene, need RRA Score < 0.05 AND best CCT FDR < 0.05
  # This requires per-gene pairing. Simplified: intersection of significant sets
  rra_genes <- unique(if(!is.null(rra) && nrow(rra)>0) rra$Name[rra$Score<ALPHA] else character(0))
  cct_genes <- gene_ids[cct_fdr < ALPHA]
  dual_sig[i] <- length(intersect(rra_genes, cct_genes))
  
  if(i %% 200 == 0 || i == 1) {
    e <- as.numeric(difftime(Sys.time(),t0,units="secs"))
    cat(sprintf(" %4d/%d | %.0fs | RRA=%.1f (%.2f%%) CCT=%.1f (%.2f%%) Dual=%d | ETA %.0fs\n",
                i, N_PERM, e, mean(rra_sig[1:i]), mean(rra_sig[1:i])/n_genes*100,
                mean(cct_sig[1:i]), mean(cct_sig[1:i])/n_genes*100, dual_sig[i], e/i*(N_PERM-i)))
  }
}

rra_fpr <- mean(rra_sig)/n_genes; cct_fpr <- mean(cct_sig)/n_genes; dual_fpr <- mean(dual_sig)/n_genes
cat(sprintf("\n=== RRA vs CCT Null (%d perms, pipeline CCT logic) ===\n", N_PERM))
cat(sprintf("RRA alone:  FPR=%.2f%% (%.1f sig)\n", rra_fpr*100, mean(rra_sig)))
cat(sprintf("CCT alone:  FPR=%.2f%% (%.1f sig)\n", cct_fpr*100, mean(cct_sig)))
cat(sprintf("RRA+CCT:    FPR=%.2f%% (%.1f sig)\n", dual_fpr*100, mean(dual_sig)))
if(mean(rra_sig)>0) cat(sprintf("Dual/RRA ratio: %.1f%%\n", mean(dual_sig)/mean(rra_sig)*100))

write_tsv(data.frame(iter=1:N_PERM, rra=rra_sig, cct=cct_sig, dual=dual_sig),
          file.path(BASE, "experiments/bombyx_enrichment/results/tables/rra_cct_null.tsv"))
cat("Saved.\n")
