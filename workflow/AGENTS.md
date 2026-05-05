# workflow/ — Core Pipeline Orchestration

## OVERVIEW

All application logic lives here. Three submodules: `rules/` (Snakemake DAG), `scripts/` (Python/R analysis), `modular_workflows/` (sub-workflow entry points). Nothing outside `workflow/` contains pipeline logic.

## STRUCTURE

```
workflow/
├── rules/                 # 14 Snakemake rule files (one per pipeline stage)
│   └── common.smk         # Shared config loader, sample parsing, helper functions
├── scripts/               # 14 Python + 5 R analysis scripts
└── modular_workflows/     # 10 sub-workflow Snakefiles for modular execution
```

## WHERE TO LOOK

| Task | Location | Notes |
|------|----------|-------|
| Understand pipeline flow | `rules/common.smk` → `Snakefile` (root) | common.smk defines global vars; Snakefile defines `rule all` targets |
| Add a new quantifier | `rules/quantification_<tool>.smk` + aggregation script in `scripts/` | See existing Kallisto/Salmon patterns |
| Add a new DEA method | `scripts/perform_quantifier_dea.R` | DESeq2 is fixed — methods.md mandates this |
| Debug a rule | Read rule's `shell:` or `script:` block → trace script inputs/outputs | Logs at `logs/<rule>/<sample>.log` |
| Add a modular entry point | `modular_workflows/` | Sub-Snakefiles that include subset of rules |

## CONVENTIONS

- **Rule files NEVER import from other rule files** — all shared state flows through `common.smk`
- **`common.smk` is the only shared dependency** — defines `SAMPLES`, `ALIGNER`, path constructors, wildcards
- **Scripts receive Snakemake context** via `snakemake` object (Python) or `snakemake@` (R) — never via CLI args
- **One Conda env per rule module** — defined as YAML in `envs/`, referenced as `conda: "../../envs/<name>.yaml"`
- **Output paths NEVER collide across rules** — `results/NN.stage/` prefix enforces this
- **Chinese comments in rule files** — rule shells and descriptions use Chinese, Snakemake DSL keywords are English

## ANTI-PATTERNS

- Don't read `config.yaml` directly in scripts — use `snakemake.config`
- Don't hardcode paths in scripts — receive them via `snakemake.input` / `snakemake.output`
- Don't create `workflow/__init__.py` — this is not a Python package
- Don't import between scripts — they run standalone per Snakemake rule invocation
