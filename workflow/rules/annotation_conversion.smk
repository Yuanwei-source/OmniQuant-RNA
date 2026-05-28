rule annotation_format_conversion:
    """
    Master rule to ensure both GTF and GFF3 formats are available
    """
    input:
        original=REFERENCE_SOURCE_ANNOTATION
    output:
        generated=(REFERENCE_GTF if REFERENCE_SOURCE_FORMAT == "gff3" else REFERENCE_ANNOTATION),
        converted = f"{REFERENCE_DIR}/annotation_conversion_complete.flag"
    run:
        input_format = REFERENCE_SOURCE_FORMAT

        if input_format == "gff3":
            shell("gffread -T {input.original} -o {output.generated}")
        else:
            shell("gffread {input.original} -o {output.generated}")

        with open(output.converted, 'w') as f:
            f.write(f"Input annotation: {input.original}\n")
            f.write(f"Detected format: {input_format}\n")
            f.write(f"Canonical GTF: {REFERENCE_GTF}\n")
            f.write(f"Canonical GFF3: {REFERENCE_ANNOTATION}\n")
            f.write(f"Transcriptome target: {REFERENCE_TRANSCRIPTOME}\n")

rule extract_transcriptome:
    """
    Extract transcriptome sequences from genome using annotation
    """
    input:
        genome=REFERENCE_GENOME,
        gtf=REFERENCE_GTF,
        conversion_flag=f"{REFERENCE_DIR}/annotation_conversion_complete.flag"
    output:
        transcriptome=REFERENCE_TRANSCRIPTOME
    conda:
        "../../envs/quantification.yaml"
    log:
        "logs/extract_transcriptome.log"
    shell:
        """
        tmp_genome=$(mktemp --suffix=.fa)
        trap 'rm -f "$tmp_genome"' EXIT
        awk '/^>/ {{sub(/ .*/, "", $0)}} {{print}}' {input.genome} > "$tmp_genome"
        gffread -w {output.transcriptome} -g "$tmp_genome" {input.gtf} > {log} 2>&1
        """
