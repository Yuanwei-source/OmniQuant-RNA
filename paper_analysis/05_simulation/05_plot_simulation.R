#!/usr/bin/env Rscript
# Figure X: Consensus performance under clean and isoform-switching stress simulations
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(tidyr); library(patchwork)
})

source("../theme_nature.R")

CLEAN <- "experiments/bombyx_enrichment/results/polyester_stress/clean_5rep/summary/summary.tsv"
STRESS <- "experiments/bombyx_enrichment/results/polyester_stress/stress_5rep/summary/summary.tsv"
STRESS_ALL <- "experiments/bombyx_enrichment/results/polyester_stress/stress_5rep/summary/all_metrics.tsv"
FIG_DIR <- "benchmark_results/figures"

clean <- read_tsv(CLEAN, show_col_types=FALSE)
stress <- read_tsv(STRESS, show_col_types=FALSE)
stress_all <- read_tsv(STRESS_ALL, show_col_types=FALSE)

methods_order <- c("Salmon","Kallisto","featureCounts","Fisher(RRA+CCT)","RRA+CCT(TierA)")
method_labels <- c("Salmon","Kallisto","FC","Fisher","Tier A")
method_fills <- c("Salmon"=col_salmon,"Kallisto"=col_kallisto,
  "featureCounts"=col_featurecounts,"Fisher(RRA+CCT)"=col_stringtie,
  "RRA+CCT(TierA)"=col_consensus)
clean$Method <- factor(clean$Method, methods_order)
stress$Method <- factor(stress$Method, methods_order)

# ── Panel A: Clean F1 ──
pA <- ggplot(clean %>% filter(!Method %in% c("CCT","RRA")), 
       aes(x=Method, y=F1_mean, fill=Method)) +
  geom_bar(stat="identity", width=0.6) +
  geom_errorbar(aes(ymin=F1_mean-F1_sd, ymax=F1_mean+F1_sd), width=0.15) +
  scale_fill_manual(values=method_fills) +
  labs(title="Clean benchmark: F1 score",y="F1",x="") +
  theme(legend.position="none") + ylim(0,1)

# ── Panel B: Stress Precision ──
pB <- ggplot(stress %>% filter(!is.na(precision_mean)), 
       aes(x=Method, y=precision_mean, fill=Method)) +
  geom_bar(stat="identity", width=0.6) +
  geom_errorbar(aes(ymin=precision_mean-precision_sd, ymax=precision_mean+precision_sd), width=0.15) +
  scale_fill_manual(values=method_fills) +
  labs(title="Stress: Precision",y="Precision",x="") +
  theme(legend.position="none") + ylim(0,0.7)

# ── Panel C: Stress Iso-switch FP% ──
pC <- ggplot(stress %>% filter(!is.na(iso_fp_rate_mean)), 
       aes(x=Method, y=iso_fp_rate_mean*100, fill=Method)) +
  geom_bar(stat="identity", width=0.6) +
  geom_errorbar(aes(ymin=(iso_fp_rate_mean-iso_fp_rate_sd)*100, 
    ymax=(iso_fp_rate_mean+iso_fp_rate_sd)*100), width=0.15) +
  scale_fill_manual(values=method_fills) +
  labs(title="Stress: Isoform-switch FP%",y="Isoform-switch false positive rate (%)",x="") +
  theme(legend.position="none")

# ── Panel D: Precision-Recall scatter ──
pD <- ggplot(stress %>% filter(!is.na(precision_mean)), 
       aes(x=recall_mean, y=precision_mean, label=Method, color=Method)) +
  geom_point(size=2.5) + geom_text(vjust=-1, size=2.5) +
  scale_color_manual(values=method_fills) +
  labs(title="Stress: Precision-Recall trade-off",x="Recall",y="Precision") +
  theme(legend.position="none") + xlim(0,0.8) + ylim(0,0.7)

# ── Panel E: Seed-paired line plot (Per-seed) ──
if ("seed" %in% names(stress_all)) {
  seed_data <- stress_all %>% 
    filter(method %in% c("Salmon","RRA+CCT(TierA)")) %>%
    select(seed, method, iso_fp_rate) %>%
    pivot_wider(names_from=method, values_from=iso_fp_rate)
  
  seed_long <- seed_data %>%
    pivot_longer(-seed, names_to="Method", values_to="iso_fp")
  seed_long$Method <- factor(seed_long$Method, c("Salmon","RRA+CCT(TierA)"))
  
  pE <- ggplot(seed_long, aes(x=Method, y=iso_fp*100, group=seed)) +
    geom_line(alpha=0.5, colour=palette_nature[["neutral_mid"]]) + 
    geom_point(size=2, colour=palette_nature[["neutral_dark"]]) +
    labs(title="Per-seed: Iso-switch FP reduction",y="Iso-switch FP%",x="") +
    scale_x_discrete(labels=c("Salmon","Tier A"))
} else {
  pE <- ggplot() + annotate("text",x=1,y=1,label="Per-seed data not available") + theme_void()
}

# ── Combine ──
design <- "
AABBCC
DDEEEE
"
combined <- pA + pB + pC + pD + pE + plot_layout(design=design, guides="collect") +
  plot_annotation(tag_levels="a") &
  theme(plot.tag=element_text(size=9, face="bold"))

save_pub_r(combined, file.path(FIG_DIR, "figure_simulation_benchmarks"),
  width_mm=183, height_mm=130)
