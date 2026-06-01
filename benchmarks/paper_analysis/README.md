# paper_analysis/ — Publication-Grade Analysis Scripts

This directory contains the FINAL, reproducible scripts used to generate
all figures and tables for the OmniQuant-RNA manuscript.

## Structure

| Directory | Content | Corresponding Figure |
|-----------|---------|---------------------|
| 01_core_benchmark/ | Annotation degradation | Fig. 2 |
| 02_decontam_p1/ | Decontamination + OFF-only | Fig. 3 |
| 03_ablation/ | Gate ablation + tier sensitivity | Fig. 4 |
| 04_stability/ | Perturbation + subsampling stability | Fig. 5 |
| 05_simulation/ | Clean + stress simulation + SEQC | Fig. 6 |
| 06_functional_coherence/ | GO/KEGG enrichment + ribosome | Fig. 7 |
| 07_figures/ | Figure assembly | — |
| 08_tables/ | Table export | — |

## Usage

```bash
./run_all.sh
```

Outputs go to `../paper_outputs/`.

## Manifest

See `00_manifest/analysis_manifest.tsv` for claim-to-script mapping.

## Original Sources

Scripts are COPIES from their original locations:
- `benchmark_results/scripts/`
- `experiments/bombyx_enrichment/scripts/`
- `experiments/drosophila_enrichment/`

Originals are preserved. These copies are the publication-frozen versions.
