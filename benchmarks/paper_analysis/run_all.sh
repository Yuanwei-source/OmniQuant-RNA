#!/usr/bin/env bash
set -euo pipefail
echo "=== OmniQuant-RNA Paper Analysis Pipeline ==="
echo ""

# All scripts run from project root (detect or use PROJECT_ROOT env var)
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/.." && pwd)}"
export PROJECT_ROOT

echo "[1/6] Core benchmark..."
python3 "$SCRIPT_DIR/01_core_benchmark/01_annotation_degradation.py"

echo "[2/6] Decontam + OFF-only..."
python3 "$SCRIPT_DIR/02_decontam_p1/01_onoff_comparison.py"

echo "[3/6] Ablation..."
python3 "$SCRIPT_DIR/03_ablation/01_gate_failure_overlap.py"
python3 "$SCRIPT_DIR/03_ablation/02_threshold_sweep.py"

echo "[4/6] Simulation..."
Rscript "$SCRIPT_DIR/05_simulation/01_clean_count_benchmark.R"
Rscript "$SCRIPT_DIR/05_simulation/02_clean_read_eval.R"
Rscript "$SCRIPT_DIR/05_simulation/03_stress_read_eval.R"

echo "[5/6] Functional coherence..."
Rscript "$SCRIPT_DIR/06_functional_coherence/01_bombyx_ora.R"
Rscript "$SCRIPT_DIR/06_functional_coherence/02_bombyx_gsea.R"
Rscript "$SCRIPT_DIR/06_functional_coherence/03_ribosome_validation.R"

echo "[6/6] Figures + tables..."
Rscript "$SCRIPT_DIR/07_figures/01_plot_functional_coherence.R"
Rscript "$SCRIPT_DIR/08_tables/export_manuscript_tables.R"

echo ""
echo "=== Done. Outputs in paper_outputs/ ==="
