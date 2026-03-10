# Quality Control Rules
# FastQC, fastp trimming, pollution removal and related QC steps

rule raw_reads_fastqc:
    input: 
        forward_reads = get_r1,
        reverse_reads = get_r2,
    output:
        forward_html = "results/01.raw_qc/{sample}_R1_fastqc.html",
        reverse_html = "results/01.raw_qc/{sample}_R2_fastqc.html",
        forward_zip = "results/01.raw_qc/{sample}_R1_fastqc.zip",
        reverse_zip = "results/01.raw_qc/{sample}_R2_fastqc.zip"
    threads: 
        config["medium"]["t"]
    resources:
        mem_mb = config["fastqc"]["memory"]
    conda:
        "../../envs/qc.yaml"
    params:
        output_dir = "results/01.raw_qc"
    log:
        "logs/{sample}_raw_fastqc.log" 
    shell:
        """
        fastqc -q -t {threads} --memory {resources.mem_mb} -o {params.output_dir} {input.forward_reads} {input.reverse_reads} >> {log} 2>&1

        forward_base=$(basename "{input.forward_reads}")
        forward_base=${{forward_base%.gz}}
        forward_base=${{forward_base%.fastq}}
        forward_base=${{forward_base%.fq}}
        reverse_base=$(basename "{input.reverse_reads}")
        reverse_base=${{reverse_base%.gz}}
        reverse_base=${{reverse_base%.fastq}}
        reverse_base=${{reverse_base%.fq}}

        mv "{params.output_dir}/${{forward_base}}_fastqc.html" "{output.forward_html}"
        mv "{params.output_dir}/${{forward_base}}_fastqc.zip" "{output.forward_zip}"
        mv "{params.output_dir}/${{reverse_base}}_fastqc.html" "{output.reverse_html}"
        mv "{params.output_dir}/${{reverse_base}}_fastqc.zip" "{output.reverse_zip}"
        """

rule quality_trimming_fastp:
    input: 
        forward_reads = get_r1,
        reverse_reads = get_r2,
        fastqc_zip = "results/01.raw_qc/{sample}_R1_fastqc.zip"
    output:
        forward_trimmed = "results/02.trimmed_data/{sample}_R1_trimmed.fastq.gz",
        reverse_trimmed = "results/02.trimmed_data/{sample}_R2_trimmed.fastq.gz",
        report_json = "results/02.trimmed_data/{sample}.fastp.json",
        report_html = "results/02.trimmed_data/{sample}.fastp.html"
    threads: 
        config["medium"]["t"]
    resources:
        mem_mb = config["medium"]["mem_mb"]
    conda:
        "../../envs/qc.yaml"
    params:
        trim_front1 = config["fastp"]["trim_front1"],
        trim_front2 = config["fastp"]["trim_front2"],
        trim_tail1 = config["fastp"]["trim_tail1"],
        trim_tail2 = config["fastp"]["trim_tail2"],
        min_length = config["fastp"]["length_required"],
        output_dir = "results/02.trimmed_data"
    log:
        "logs/{sample}_fastp_trimming.log" 
    shell:
        """
        fastp  -i {input.forward_reads} -I {input.reverse_reads} \
               -o {output.forward_trimmed} -O {output.reverse_trimmed} \
               -f {params.trim_front1} -F {params.trim_front2} \
               -t {params.trim_tail1} -T {params.trim_tail2} \
               -w {threads} \
               --detect_adapter_for_pe \
               --trim_poly_x \
               --cut_front \
               --cut_tail \
               --length_required {params.min_length} \
               -h {output.report_html} \
               -j {output.report_json} >> {log} 2>&1\
        """ 

rule trimmed_reads_fastqc:
    input: 
        forward_trimmed = "results/02.trimmed_data/{sample}_R1_trimmed.fastq.gz",
        reverse_trimmed = "results/02.trimmed_data/{sample}_R2_trimmed.fastq.gz"
    output:
        forward_html = "results/02.trimmed_data/{sample}_R1_trimmed_fastqc.html",
        reverse_html = "results/02.trimmed_data/{sample}_R2_trimmed_fastqc.html",
        forward_zip = "results/02.trimmed_data/{sample}_R1_trimmed_fastqc.zip",
        reverse_zip = "results/02.trimmed_data/{sample}_R2_trimmed_fastqc.zip"
    threads: 
        config["medium"]["t"]
    resources:
        mem_mb = config["fastqc"]["memory"]
    conda:
        "../../envs/qc.yaml"
    params:
        output_dir = "results/02.trimmed_data"
    log:
        "logs/{sample}_trimmed_fastqc.log" 
    shell:
        """
        fastqc -q -t {threads} --memory {resources.mem_mb} -o {params.output_dir} {input.forward_trimmed} {input.reverse_trimmed} >> {log} 2>&1
        """