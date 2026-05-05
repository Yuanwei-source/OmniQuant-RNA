# PROJECT KNOWLEDGE BASE

**Generated:** 2026-05-05 12:47 UTC
**Commit:** 056fe1e
**Branch:** discussion/consensus-downstream

## OVERVIEW

OmniQuant-RNA is a **Snakemake workflow** for RNA-seq quantification with multi-quantifier consensus DEA. Not a Python package — all logic lives in `workflow/`. Pipelines orchestrate bioinformatics tools via Conda environments. Host-focused analysis with optional microbial clue sidecar. Chinese primary documentation language.

## STRUCTURE

```
./
├── Snakefile              # Main pipeline DAG (includes all workflow/rules/*.smk)
├── run_analysis.sh         # Full pipeline entry script
├── run_modular.sh          # Module-based entry script (qc|alignment|quantification|dea|dry-run)
├── config/
│   └── config.yaml         # Single-source config (samples, refs, tools, DEA, consensus, decontam)
├── envs/                   # 9 per-module Conda YAML environments (.yaml)
├── workflow/               # ALL code — see workflow/AGENTS.md
│   ├── rules/              # 14 Snakemake rule modules
│   ├── scripts/            # 14 Python + 5 R per-rule scripts
│   └── modular_workflows/  # 10 sub-workflow entry points
├── docs/                   # Formal docs: methods.md, outputs.md, repository-guidelines.md
├── data/                   # Input data (FASTQ, reference) — NOT committed
├── results/                # Pipeline outputs organized by stage (00-08) — NOT committed
└── logs/                   # Per-rule Snakemake logs — NOT committed
```

## WHERE TO LOOK

| Task | Location | Notes |
|------|----------|-------|
| Add a new pipeline stage | `workflow/rules/` + register in `Snakefile` | Follow rule conventions in `workflow/rules/AGENTS.md` |
| Add a new analysis script | `workflow/scripts/` | Python for aggregation/namespace, R for DEA/plotting |
| Change tool parameters | `config/config.yaml` | Central config; scripts read via `snakemake.config` |
| Change conda dependencies | `envs/<module>.yaml` | One env per pipeline module |
| Add a modular sub-workflow | `workflow/modular_workflows/` | Create `<module>_only.smk` |
| Understand the pipeline DAG | `Snakefile` → `workflow/rules/common.smk` | common.smk defines shared vars and helper functions |
| Understand DEA/consensus | `workflow/scripts/perform_quantifier_dea.R` + `run_consensus_dea.R` | R scripts, not Python |
| Understand decontam boundary | `docs/methods.md` | Current design principles |

## CONVENTIONS

### General
- **Documentation language: Chinese** (README, comments, error messages, commit messages)
- **No `pip install`** — this is a Snakemake workflow, clone-and-run
- **No `setup.py` / `pyproject.toml`** — dependencies via `envs/*.yaml` Conda only
- **Output directories: `results/NN.stage/`** — sequential prefix enforces pipeline order
- **No tests, no CI, no Docker** — validated via `snakemake --dry-run` and manual execution
- **Dual remote**: GitHub (`origin`) + Gitee (`gitee`)

### Snakemake (see `workflow/rules/AGENTS.md`)
- `UPPER_SNAKE_CASE` for global variables defined in `common.smk`
- Rules per module: `quantification_kallisto.smk`, `alignment.smk`, etc.
- `conda: "../../envs/<module>.yaml"` per rule
- Logs at `logs/<rule>/<sample>.log`
- Shell blocks: f-string paths, Snakemake wildcards via `{wildcards.key}`

### Python (see `workflow/scripts/AGENTS.md`)
- `#!/usr/bin/env python3` shebang on standalone scripts
- `argparse` for CLI (no click/typer)
- `polars` (as `pl`) preferred over `pandas` for data processing
- `snake_case` for variables/functions; `UPPER_SNAKE_CASE` for constants
- `if __name__ == "__main__": main()` pattern
- `print()` for logging (no `logging` module)

### R (see `workflow/scripts/AGENTS.md`)
- `box::use()` for explicit package imports
- Snakemake objects: `snakemake@input`, `snakemake@output`, `snakemake@params`, `snakemake@config`
- `sink(snakemake@log)` for log redirection
- `dplyr` + `readr` for data; `DESeq2` for DEA (fixed — no edgeR/limma)

### Config
- Single `config/config.yaml` — no multi-file configs
- Comparisons: `"treatment_vs_control"` format (left=treatment, right=control)
- `config["key"]` with `dict.get("key", default)` pattern in scripts

### Docs
- Formal docs: `README.md` + 3 `docs/*.md` only (4 files total)
- Temporary/audit/draft content → `experiments/` (local only, not committed)
- Do NOT add new formal docs without changing the 3-mandated-to-4 balance

## ANTI-PATTERNS (THIS PROJECT)

- **NOT a Python package** — don't add `setup.py`, `pyproject.toml`, or `__init__.py`
- **Don't add `requirements.txt`** — use Conda YAML in `envs/` only
- **Don't add edgeR or limma** — DEA is DESeq2-only
- **Don't change consensus stats** — RRA + CCT as dual co-primary engines; no dynamic weighting
- **Don't mix host and microbe results** in single DEG tables — they're kept separate
- **Don't commit `data/`, `results/`, `logs/`, `benchmarks/`** — .gitignore excludes them
- **Don't create new formal docs** unless existing 3 can't absorb the content
- **Don't use `logging` module in Python scripts** — use `print()`
- **Don't assume any linter/formatter config** — none exist (no black, ruff, flake8, pre-commit)

## UNIQUE STYLES

- **English Snakemake rules + Chinese comments/docs** — mixed language is normal here
- **`polars` over `pandas`** — this project specifically prefers polars for data processing
- **Per-rule Conda envs** — not one monolithic env; isolation per pipeline stage
- **Decontam has its own clue sidecar** — `results/03.decontam/clues/` is SEPARATE from the host DEA pipeline
- **Reference auto-detection** — `reference_config.smk` scans `data/reference/` for FASTA/annotation files
- **No versioning** — no git tags, no `version.py`, no changelog

## COMMANDS

```bash
# Full pipeline
./run_analysis.sh
CORES=32 nohup ./run_analysis.sh > analysis.log 2>&1 &

# Modular execution
./run_modular.sh qc --cores 8
./run_modular.sh quantification --cores 16
./run_modular.sh dea --cores 16
./run_modular.sh all --cores 24

# Direct Snakemake
snakemake --use-conda --cores 32 --printshellcmds
snakemake -s workflow/modular_workflows/quantification_only.smk --use-conda --cores 16

# Validation
snakemake --dry-run
./run_modular.sh dry-run
snakemake --rulegraph | dot -Tpng > rules_graph.png

# Sample prep
python workflow/scripts/generate_samples.py data/fastq/ -o data/fastq/samples.tsv
```

## NOTES

- `envs/dea.yaml` uses default conda channels while others use Tsinghua mirrors — be aware of resolution differences
- Some env YAML files have duplicate dependency pins (e.g., `xopen>=1.7` AND `xopen>=2.0` in `decontam.yaml`)
- `data/reference/` must contain `genome.fasta` + annotation — autodetected by `reference_config.smk`
- Branch `discussion/consensus-downstream` suggests active consensus DEA discussion
- `run_analysis.sh` self-bootstraps conda if `snakemake` not on PATH; expects a `snakemake` conda env
