# workflow/modular_workflows/ — Sub-Workflow Entry Points

## OVERVIEW

10 standalone Snakefiles that dispatch subsets of the full pipeline for modular execution. Each includes only the rule files needed for its target stage. Invoked via `./run_modular.sh <module>` or directly with `snakemake -s`.

## FILES

| File | Module | Includes |
|------|--------|----------|
| `qc_only.smk` | QC + decontam | qc.smk, decontam.smk |
| `decontam_only.smk` | Decontam only | decontam.smk |
| `alignment_only.smk` | Alignment only | common.smk, alignment.smk, decontam.smk |
| `quantification_only.smk` | All quantifiers | common.smk, 4 quantification rules |
| `quantification_kallisto_only.smk` | Kallisto only | common.smk, quantification_kallisto.smk |
| `quantification_salmon_only.smk` | Salmon only | common.smk, quantification_salmon.smk |
| `quantification_featurecounts_only.smk` | featureCounts only | common.smk, quantification_featurecounts.smk |
| `quantification_stringtie_only.smk` | StringTie only | common.smk, quantification_stringtie.smk |
| `differential_expression_only.smk` | DEA only | common.smk, differential_expression.smk, consensus_expression.smk |
| `annotation_conversion_only.smk` | Annotation prep | annotation_conversion.smk |

## PATTERN

Each file follows the same template:
```python
configfile: "../../config/config.yaml"
include: "../../workflow/rules/common.smk"
include: "../../workflow/rules/<required_rules>.smk"

rule all:
    input:
        # List of final output files for this sub-workflow
```

## HOW THEY MAP TO run_modular.sh

```
./run_modular.sh qc              → qc_only.smk
./run_modular.sh decontam        → decontam_only.smk
./run_modular.sh alignment       → alignment_only.smk
./run_modular.sh quantification  → quantification_only.smk
./run_modular.sh dea             → differential_expression_only.smk
./run_modular.sh all             → root Snakefile
```

Single-quantifier targets are available via direct Snakemake invocation:
```bash
snakemake -s workflow/modular_workflows/quantification_kallisto_only.smk --use-conda --cores 16
```

## ANTI-PATTERNS

- Don't add new modular files for stages that already exist — extend the existing sub-workflow
- Don't define rules in modular files — only `include:` existing rule files and `rule all:` targets
- Don't modify `configfile:` path — always `"../../config/config.yaml"`
- Don't skip `common.smk` — all modular files depend on its variable definitions
