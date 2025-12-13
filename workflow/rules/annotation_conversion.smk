# Annotation Conversion Rules
# Convert between GFF3 and GTF formats

import os

def detect_annotation_format(annotation_file):
    """
    Detect whether the annotation file is GFF3 or GTF format
    """
    if annotation_file.endswith('.gff3') or annotation_file.endswith('.gff'):
        return 'gff3'
    elif annotation_file.endswith('.gtf'):
        return 'gtf'
    else:
        # Read first few lines to detect format
        try:
            with open(annotation_file, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    if line.strip() == '':
                        continue
                    fields = line.split('\t')
                    if len(fields) >= 9:
                        # Check 9th column format
                        if 'gene_id' in fields[8] and 'transcript_id' in fields[8]:
                            return 'gtf'
                        elif '=' in fields[8]:
                            return 'gff3'
                    break
        except:
            pass
        # Default to gff3 if can't detect
        return 'gff3'

# Determine input format and set conversion target
input_annotation = config["reference"]["gff3"]
input_format = detect_annotation_format(input_annotation)

if input_format == 'gff3':
    target_format = 'gtf'
    target_file = input_annotation.replace('.gff3', '.gtf').replace('.gff', '.gtf')
else:
    target_format = 'gff3'
    target_file = input_annotation.replace('.gtf', '.gff3')

rule convert_gff3_to_gtf:
    input:
        gff3 = config["reference"]["gff3"]
    output:
        gtf = "data/reference/genome.gtf"
    conda:
        "../../envs/quantification.yaml"
    log:
        "logs/gff3_to_gtf_conversion.log"
    shell:
        """
        # Use gffread to convert GFF3 to GTF
        gffread -T {input.gff3} -o {output.gtf} >> {log} 2>&1
        """

rule convert_gtf_to_gff3:
    input:
        gtf = lambda wildcards: config["reference"]["gff3"] if config["reference"]["gff3"].endswith('.gtf') else "data/reference/genome.gtf"
    output:
        gff3 = "data/reference/genome_converted.gff3"
    conda:
        "../../envs/quantification.yaml"
    log:
        "logs/gtf_to_gff3_conversion.log"
    shell:
        """
        # Use gffread to convert GTF to GFF3
        gffread {input.gtf} -o {output.gff3} >> {log} 2>&1
        """

rule create_annotation_symlinks:
    """
    Create symbolic links for both GTF and GFF3 formats
    This ensures both formats are available regardless of input
    """
    input:
        original = config["reference"]["gff3"]
    output:
        gtf_link = "data/reference/annotation.gtf",
        gff3_link = "data/reference/annotation.gff3"
    run:
        import os
        
        # Detect input format
        input_format = detect_annotation_format(input.original)
        
        if input_format == 'gff3':
            # Create GTF by conversion, link GFF3 directly
            shell("gffread -T {input.original} -o data/reference/genome.gtf")
            os.symlink(os.path.basename(input.original), output.gff3_link)
            os.symlink("genome.gtf", output.gtf_link)
        else:
            # Create GFF3 by conversion, link GTF directly
            shell("gffread {input.original} -o data/reference/genome.gff3")
            os.symlink(os.path.basename(input.original), output.gtf_link)
            os.symlink("genome.gff3", output.gff3_link)

rule annotation_format_conversion:
    """
    Master rule to ensure both GTF and GFF3 formats are available
    """
    input:
        original = config["reference"]["gff3"]
    output:
        converted = "data/reference/annotation_conversion_complete.flag"
    run:
        import os
        
        # Detect input format
        input_format = detect_annotation_format(input.original)
        
        if input_format == 'gff3':
            # Convert to GTF
            output_file = "data/reference/genome.gtf"
            shell("gffread -T {input.original} -o {output_file}")
        else:
            # Convert to GFF3  
            output_file = "data/reference/genome.gff3"
            shell("gffread {input.original} -o {output_file}")
        
        # Create completion flag
        with open(output.converted, 'w') as f:
            f.write(f"Conversion completed: {input_format} -> {'gtf' if input_format == 'gff3' else 'gff3'}\n")
            f.write(f"Input file: {input.original}\n")
            f.write(f"Output file: {output_file}\n")

rule extract_transcriptome:
    """
    Extract transcriptome sequences from genome using annotation
    """
    input:
        genome = config["reference"]["genome"],
        gtf = "data/reference/genome.gtf",
        conversion_flag = "data/reference/annotation_conversion_complete.flag"
    output:
        transcriptome = config["reference"]["transcriptome"]
    conda:
        "../../envs/quantification.yaml"
    log:
        "logs/extract_transcriptome.log"
    shell:
        """
        gffread -w {output.transcriptome} -g {input.genome} {input.gtf} 2> {log}
        """
