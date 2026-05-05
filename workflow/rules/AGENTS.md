# workflow/rules/ — Snakemake Rule Definitions

## OVERVIEW

14 Snakemake rule files, each defining one pipeline stage. Included by root `Snakefile` via `include:`. `common.smk` is the only shared dependency — it loads config, parses samples, and defines global variables.

## FILES

| File | Stage | Key Output |
|------|-------|------------|
| `common.smk` | Config loader + shared helpers | `SAMPLES`, `ALIGNER`, path constructors, decontam targets |
| `reference_config.smk` | Reference auto-detection | Symlinks `genome.fasta`, `genome.gff3` |
| `annotation_conversion.smk` | GFF3 normalization | `featurecounts.gtf` |
| `qc.smk` | FastQC + fastp | Trimmed FASTQ |
| `decontam.smk` | Host removal + Kraken2 | Clean FASTQ + audit + microbial clues |
| `alignment.smk` | HISAT2 or STAR | BAM + BAI |
| `reference_namespace.smk` | Unified tx2gene / gene namespace | `tx2gene_master.tsv`, `gene_namespace.tsv`, import manifests |
| `quantification_kallisto.smk` | Kallisto quant + matrix | `abundance.tsv` → `kallisto_*_matrix.tsv` |
| `quantification_salmon.smk` | Salmon quant + matrix | `quant.sf` → `salmon_*_matrix.tsv` |
| `quantification_featurecounts.smk` | featureCounts quant + matrix | `counts.txt` → `featurecounts_*_matrix.tsv` |
| `quantification_stringtie.smk` | StringTie quant + merge + matrix | `gene_abundances.tab` → `stringtie_*_matrix.tsv` |
| `differential_expression.smk` | DESeq2 per quantifier | `results/06.differential_expression/<quant>/` |
| `consensus_expression.smk` | RRA + CCT consensus | `results/07.consensus_expression/<contrast>/` |
| `report.smk` | MultiQC | `results/08.reports/multiqc_report.html` |

## RULE CONVENTIONS

```
rule <name>:
    input: <files or function returning list>
    output: <files>
    conda: "../../envs/<module>.yaml"
    log: "logs/<rule>/<sample>.log"
    threads: <config value>
    shell: """ ... """  # or script: "..."
```

- **Global vars**: Defined in `common.smk`, referenced as `{VAR_NAME}` in rules — e.g. `KALLISTO_NATIVE_DIR`, `SAMPLES`
- **Path construction**: Use f-strings or `.format()` with `{wildcards.key}` and `{VAR_NAME}`
- **Input functions**: `def get_input(wildcards):` pattern when input depends on wildcards
- **`ruleorder:`** used to resolve ambiguous output paths (e.g., `hisat2_align > select_alignment_bam`)
- **`checkpoint:`** not used — all rules are deterministic
- **Conditional targets**: `get_decontam_all_targets(SAMPLES)` in `common.smk` returns empty list when decontam disabled
- **Named target rules**: `rule qc_only:`, `rule stringtie_only:`, etc. for modular execution
- **`run:` blocks**: Rarely used — prefer `shell:` or `script:`

## ANTI-PATTERNS

- Don't put logic in rule files beyond input/output/threads/conda — logic goes in scripts
- Don't use `module:` (Snakemake v8+) — the codebase uses `include:`
- Don't hardcode paths in `output:` — use `get_path_constructor()` from `common.smk`
- Don't duplicate `input:` definitions — share them via `common.smk` helpers
