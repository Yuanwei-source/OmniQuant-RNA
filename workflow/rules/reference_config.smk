import os

REFERENCE_DIR = "data/reference"
DEFAULT_REFERENCE_GENOME = os.path.join(REFERENCE_DIR, "genome.fasta")
DEFAULT_REFERENCE_ANNOTATION = os.path.join(REFERENCE_DIR, "genome.gff3")
DEFAULT_REFERENCE_GTF = os.path.join(REFERENCE_DIR, "genome.gtf")
DEFAULT_REFERENCE_GFF3 = os.path.join(REFERENCE_DIR, "genome.gff3")
DEFAULT_REFERENCE_TRANSCRIPTOME = os.path.join(REFERENCE_DIR, "transcriptome.fasta")
LEGACY_REFERENCE_KEYS = ("gff3", "gtf", "transcriptome")

REFERENCE_AUTODETECT_ANNOTATION_PLACEHOLDERS = {
    None,
    "",
    DEFAULT_REFERENCE_ANNOTATION,
    DEFAULT_REFERENCE_GTF,
    DEFAULT_REFERENCE_GFF3,
    os.path.join(REFERENCE_DIR, "annotation.gtf"),
    os.path.join(REFERENCE_DIR, "annotation.gff3"),
}


def detect_annotation_format(annotation_file):
    lower_file = str(annotation_file).lower()
    if lower_file.endswith((".gtf", ".gtf.gz")):
        return "gtf"
    if lower_file.endswith((".gff", ".gff3", ".gff.gz", ".gff3.gz")):
        return "gff3"
    raise ValueError(
        f"Unsupported annotation extension for '{annotation_file}'. Use .gff, .gff3, or .gtf"
    )


def maybe_symlink_reference(source_path, target_path):
    source_abs = os.path.abspath(source_path)
    target_abs = os.path.abspath(target_path)

    if source_abs == target_abs:
        return False

    os.makedirs(os.path.dirname(target_abs), exist_ok=True)

    if os.path.lexists(target_abs):
        if os.path.realpath(target_abs) == source_abs:
            return False
        if os.path.islink(target_abs) or os.path.isfile(target_abs):
            os.remove(target_abs)
        else:
            raise ValueError(f"Refusing to replace non-file path: {target_abs}")

    os.symlink(source_abs, target_abs)
    return True


def auto_detect_references():
    reference_config = config.setdefault("reference", {})
    if not os.path.exists(REFERENCE_DIR):
        return

    if reference_config.get("genome", DEFAULT_REFERENCE_GENOME) == DEFAULT_REFERENCE_GENOME:
        fasta_files = [
            file_name for file_name in os.listdir(REFERENCE_DIR)
            if file_name.endswith((".fa", ".fasta", ".fna"))
            and file_name not in ["genome.fasta", "transcriptome.fasta"]
            and not file_name.endswith(".tmp")
            and not file_name.endswith(".clean")
        ]
        if len(fasta_files) == 1:
            source_path = os.path.join(REFERENCE_DIR, fasta_files[0])
            if maybe_symlink_reference(source_path, DEFAULT_REFERENCE_GENOME):
                print(
                    f"\n[OmniQuant-RNA] Auto-detected: symlinked '{fasta_files[0]}' to 'genome.fasta'"
                )
            reference_config["genome"] = DEFAULT_REFERENCE_GENOME

    annotation_value = reference_config.get("annotation")
    if annotation_value in REFERENCE_AUTODETECT_ANNOTATION_PLACEHOLDERS:
        annotation_files = [
            file_name for file_name in os.listdir(REFERENCE_DIR)
            if file_name.endswith((".gff", ".gff3", ".gtf"))
            and file_name not in [
                "genome.gff3",
                "genome.gtf",
                "annotation.gtf",
                "annotation.gff3",
                "genome.featurecounts.gtf",
            ]
            and not file_name.startswith("genome_converted")
        ]
        if len(annotation_files) == 1:
            source_path = os.path.join(REFERENCE_DIR, annotation_files[0])
            reference_config["annotation"] = source_path


auto_detect_references()

REFERENCE_CONFIG = config.setdefault("reference", {})
LEGACY_REFERENCE_KEYS_PRESENT = [key for key in LEGACY_REFERENCE_KEYS if key in REFERENCE_CONFIG]
REFERENCE_GENOME = REFERENCE_CONFIG.get("genome", DEFAULT_REFERENCE_GENOME)
REFERENCE_SOURCE_ANNOTATION = REFERENCE_CONFIG.get("annotation")

if not REFERENCE_SOURCE_ANNOTATION:
    if LEGACY_REFERENCE_KEYS_PRESENT:
        raise KeyError(
            "Legacy reference keys gff3/gtf/transcriptome are no longer supported. "
            "Use reference.annotation and let the workflow derive GTF/GFF3/transcriptome internally."
        )
    raise KeyError(
        "Missing reference annotation. Set reference.annotation to a .gff, .gff3, or .gtf file."
    )

REFERENCE_SOURCE_FORMAT = detect_annotation_format(REFERENCE_SOURCE_ANNOTATION)
REFERENCE_GTF = (
    REFERENCE_SOURCE_ANNOTATION if REFERENCE_SOURCE_FORMAT == "gtf" else DEFAULT_REFERENCE_GTF
)
REFERENCE_GFF3 = (
    REFERENCE_SOURCE_ANNOTATION if REFERENCE_SOURCE_FORMAT == "gff3" else DEFAULT_REFERENCE_GFF3
)
REFERENCE_TRANSCRIPTOME = DEFAULT_REFERENCE_TRANSCRIPTOME

REFERENCE_CONFIG["annotation"] = REFERENCE_SOURCE_ANNOTATION
