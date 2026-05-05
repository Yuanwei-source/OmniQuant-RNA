# workflow/scripts/ — Analysis Scripts

## OVERVIEW

19 standalone scripts (14 Python + 5 R) invoked by Snakemake rules. Python handles data aggregation, namespace construction, and manifest generation. R handles all statistical analysis (DESeq2 DEA, RRA/CCT consensus, visualization).

## LANGUAGE SPLIT

| Language | Count | Purpose |
|----------|-------|---------|
| **Python** | 14 | Aggregation (featureCounts, Salmon, Kallisto, StringTie), namespace building, manifest generation, microbe clue summary |
| **R** | 5 | DESeq2 DEA (`perform_quantifier_dea.R`), consensus aggregation (`run_consensus_dea.R`), result integration, microbe plotting, host context overlay |

## PYTHON CONVENTIONS

### Boilerplate
```python
#!/usr/bin/env python3
import argparse
import polars as pl  # preferred over pandas

def main():
    # parse args, do work, write output

if __name__ == "__main__":
    main()
```

### Key Rules
- **Library priority**: `polars` (as `pl`) > `pandas` (as `pd`) — use polars unless blocked
- **CLI**: `argparse` only — no click, typer, or fire
- **Snakemake access**: `snakemake.input`, `snakemake.output`, `snakemake.params`, `snakemake.config`
- **Logging**: `print()` — no `logging` module
- **Paths**: `pathlib.Path` preferred
- **Shebang**: `#!/usr/bin/env python3` on ALL standalone scripts
- **Error handling**: `raise ValueError(...)` for data issues, `sys.exit(1)` for CLI errors
- **Output format**: TSV (tab-separated) — `.tsv` extension
- **Types**: Optional type annotations on shared functions; scripts skip them

### Script Categories
- **Aggregation** (`aggregate_*.py`): Read per-sample native quantifier outputs → produce unified count/TPM matrices
- **Namespace** (`build_*.py`): Construct tx2gene, gene namespace, StringTie bridge, import manifests
- **Utilities** (`generate_samples.py`, `normalize_featurecounts_gtf.py`, `create_gene_mapping.py`, `summarize_microbe_clues.py`)

## R CONVENTIONS

### Boilerplate
```r
box::use(dplyr[...], readr[...], ...)
sink(snakemake@log)

# Read inputs from snakemake@input
data <- readr::read_tsv(snakemake@input[["matrix"]])

# ... analysis ...

# Write outputs to snakemake@output
readr::write_tsv(result, snakemake@output[["results"]])
sink()
```

### Key Rules
- **Package imports**: `box::use()` for explicit namespace — `dplyr`, `readr`, `DESeq2`, `ggplot2`, `tximport`
- **Snakemake access**: `snakemake@input`, `snakemake@output`, `snakemake@params`, `snakemake@config`, `snakemake@log`, `snakemake@threads`
- **Logging**: `sink(snakemake@log)` at top — redirects all output
- **DEA**: `DESeq2` only — no edgeR, no limma-voom
- **Consensus**: `RobustRankAggreg` for RRA, custom CCT (Cauchy Combination Test) — no dynamic weighting
- **Lint suppression**: `# nolint start: object_usage_linter.` / `# nolint end` for known false positives
- **Helper operators**: `%||%` (null coalesce), `as_flag()` (boolean parse)

## ANTI-PATTERNS

- Don't add `__init__.py` — scripts are standalone, not a package
- Don't import between scripts — each runs independently per Snakemake rule
- Don't add click/typer/fire — argparse only for Python
- Don't use `edgeR` or `limma` in R — DEA is DESeq2-fixed
- Don't add new consensus stat methods without updating `docs/methods.md`
- Don't use `logging` module — `print()` only
