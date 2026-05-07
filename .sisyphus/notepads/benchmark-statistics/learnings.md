# Benchmark Statistics Script — Learnings

## Data Structures
- `annotation_degradation.tsv`: mode, level, method, global_recall, seed — 5 seeds per (mode×level) combination
  - salmon/kallisto have constant global_recall=1.0 at all levels/modes (annotation-independent quantifiers)
  - At 0.5 level, many modes show constant values across seeds (single unique value), making paired t-test impossible
- `subsampling_stability.tsv`: gene-level stability fractions (0-1) per method; ~4000 genes, many with 0 stability
- `ablation_summary.tsv`: 4 config rows (baseline, namespace, consensus, decontam) with summary metrics

## Constant Data Handling
- When paired differences have sd ≈ 0, `t.test()` fails with "data are essentially constant"
- Need to check `sd(c_vals - q_vals) < 1e-10` before calling t.test
- For constant data, report "— (constant diff)" as effect size and "—" for CI/p-value

## R Conventions in This Project
- Shebang: `#!/usr/bin/env Rscript`
- Packages: `suppressPackageStartupMessages({library(readr); library(dplyr); library(tibble)})`
- File reading: `read_tsv(..., show_col_types = FALSE)`
- Error handling: use `tryCatch()` around t.test, wilcox.test
- Output: TSV via `write_tsv()`
