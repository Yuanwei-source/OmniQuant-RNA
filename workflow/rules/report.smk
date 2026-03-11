# Report Generation Rules
# MultiQC and custom reports

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
        # Quantification results
        kallisto=expand("results/04.quantification/native/kallisto/per_sample/{sample}/run_info.json", sample=SAMPLES)
    output:
        "results/07.reports/multiqc_report.html"
    conda:
        "../../envs/qc.yaml"
    log:
        "logs/multiqc.log"
    shell:
        """
        multiqc results/01.raw_qc/ results/02.trimmed_data/ results/04.quantification/native/kallisto/ results/04.quantification/native/salmon/ \
        -o results/07.reports/ -f \
        --filename multiqc_report.html \
        --title "RNA-seq Comprehensive Quality Control Report" \
        --comment "Complete analysis from raw reads through quantification" 2> {log}
        """
