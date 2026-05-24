import os
import sys

# =============================================================================
# OmniQuant-RNA Reference Configuration
# 
# All reference paths are derived from a single config key: reference.species
# 
# Required directory structure:
#   data/reference/{species}/
#   ├── genome.fasta          ← MUST exist (symlink ok)
#   ├── annotation.gff3       ← MUST exist (symlink ok)
#   ├── genome.gtf            ← pipeline-generated
#   ├── genome.featurecounts.gtf  ← pipeline-generated
#   ├── transcriptome.fasta   ← pipeline-generated
#   ├── annotation_conversion_complete.flag  ← pipeline-generated
#   ├── hisat2_index/         ← pipeline-generated
#   ├── salmon_index/         ← pipeline-generated
#   └── kallisto_index/       ← pipeline-generated
# =============================================================================

# ── Species → Reference Directory ────────────────────────────────────────────
REFERENCE_SPECIES = config["reference"]["species"]
REFERENCE_DIR = f"data/reference/{REFERENCE_SPECIES}"

# ── Input files (MUST exist) ─────────────────────────────────────────────────
REFERENCE_GENOME = f"{REFERENCE_DIR}/genome.fasta"
REFERENCE_ANNOTATION = f"{REFERENCE_DIR}/annotation.gff3"

if not os.path.isfile(REFERENCE_GENOME):
    raise ValueError(f"Missing genome fasta: {REFERENCE_GENOME}")
if not os.path.isfile(REFERENCE_ANNOTATION):
    raise ValueError(f"Missing annotation gff3: {REFERENCE_ANNOTATION}")

# ── Derived files (pipeline-generated) ───────────────────────────────────────
REFERENCE_GTF = f"{REFERENCE_DIR}/genome.gtf"
FEATURECOUNTS_GTF = f"{REFERENCE_DIR}/genome.featurecounts.gtf"
TRANSCRIPTOME_FASTA = f"{REFERENCE_DIR}/transcriptome.fasta"

# ── Index directories (pipeline-generated) ───────────────────────────────────
HISAT2_INDEX_DIR = f"{REFERENCE_DIR}/hisat2_index"
SALMON_INDEX_DIR = f"{REFERENCE_DIR}/salmon_index"
KALLISTO_INDEX_DIR = f"{REFERENCE_DIR}/kallisto_index"

# ── Backward-compatible aliases ──────────────────────────────────────────────
REFERENCE_SOURCE_ANNOTATION = REFERENCE_ANNOTATION
REFERENCE_TRANSCRIPTOME = TRANSCRIPTOME_FASTA

# ── Annotation format detection ──────────────────────────────────────────────
if REFERENCE_ANNOTATION.endswith(".gtf"):
    REFERENCE_GTF = REFERENCE_ANNOTATION  # source IS gtf, no conversion needed
    REFERENCE_SOURCE_FORMAT = "gtf"
else:
    REFERENCE_SOURCE_FORMAT = "gff3"
