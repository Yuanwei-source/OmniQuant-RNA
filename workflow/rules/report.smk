# Report Generation Rules
# MultiQC and custom reports


def get_multiqc_optional_inputs(samples):
    if not decontam_enabled():
        return []

    return [
        expand("results/03.decontam/qc/{sample}_{read}_clean_fastqc.zip", sample=samples, read=["R1", "R2"]),
        "results/03.decontam/stats/project_decontam_summary.tsv",
    ]


def get_multiqc_scan_targets():
    scan_targets = [
        "results/01.raw_qc/",
        "results/02.trimmed_data/",
        "results/05.quantification/native/kallisto/",
        "results/05.quantification/native/salmon/",
    ]
    if decontam_enabled():
        scan_targets.append("results/03.decontam/")
    return " ".join(scan_targets)

rule multiqc:
    """
    Generate comprehensive quality control report with MultiQC
    Including raw QC, trimmed data, clean data, and quantification results
    """
    input:
        # Raw QC results
        raw_fastqc=expand("results/01.raw_qc/{sample}_{read}_fastqc.zip", sample=SAMPLES, read=["R1", "R2"]),
        # Trimmed data QC results  
        fastp_reports=expand("results/02.trimmed_data/{sample}.fastp.json", sample=SAMPLES),
        trimmed_fastqc=expand("results/02.trimmed_data/{sample}_{read}_trimmed_fastqc.zip", sample=SAMPLES, read=["R1", "R2"]),
        decontam_optional=get_multiqc_optional_inputs(SAMPLES),
        # Quantification results
        kallisto=expand("results/05.quantification/native/kallisto/per_sample/{sample}/run_info.json", sample=SAMPLES)
    output:
        "results/08.reports/multiqc_report.html"
    conda:
        "../../envs/qc.yaml"
    params:
        scan_targets=get_multiqc_scan_targets()
    log:
        "logs/multiqc.log"
    shell:
        """
        multiqc {params.scan_targets} \
        -o results/08.reports/ -f \
        --filename multiqc_report.html \
        --title "RNA-seq Comprehensive Quality Control Report" \
        --comment "Complete analysis from raw reads through quantification" 2> {log}
        """
