DECONTAM_DIR = "results/02.5.decontam"
DECONTAM_TMP_DIR = f"{DECONTAM_DIR}/tmp"
DECONTAM_CLEAN_DIR = f"{DECONTAM_DIR}/clean"
DECONTAM_AUDIT_DIR = f"{DECONTAM_DIR}/audit"
DECONTAM_STATS_DIR = f"{DECONTAM_DIR}/stats"
DECONTAM_QC_DIR = f"{DECONTAM_DIR}/qc"
DECONTAM_REF_DIR = f"{DECONTAM_DIR}/reference"
DECONTAM_INDEX_DIR = f"{DECONTAM_DIR}/index"
DECONTAM_CONFIG = config.get("decontam", {})


def decontam_audit_enabled():
    return bool(DECONTAM_CONFIG.get("outputs", {}).get("keep_audit_fastq", True))


def get_decontam_reference(key, default=""):
    value = DECONTAM_CONFIG.get("references", {}).get(key, default)
    if value in (None, "null"):
        return ""
    return value


HOST_RESCUE_EXTRA = DECONTAM_CONFIG.get("host_rescue", {}).get("extra", "")
HOST_RESCUE_MAX_INSERT = DECONTAM_CONFIG.get("host_rescue", {}).get("max_insert", 1000)
PRESCREEN_EXTRA = DECONTAM_CONFIG.get("prescreen", {}).get("extra", "")
PRESCREEN_MAX_INSERT = DECONTAM_CONFIG.get("prescreen", {}).get("max_insert", 1000)
CLASSIFIER_CONFIG = DECONTAM_CONFIG.get("classifier", {})
CLASSIFIER_DB = CLASSIFIER_CONFIG.get("db", "")
CLASSIFIER_THREADS = max(1, int(CLASSIFIER_CONFIG.get("threads", config.get("medium", {}).get("t", 4))))


def get_reference_list(key):
    value = DECONTAM_CONFIG.get("references", {}).get(key, [])
    if value in (None, "null", ""):
        return []
    if isinstance(value, str):
        return [value]
    return list(value)


rule decontam_audit_reference_genome:
    """
    Pre-check the host reference genome to identify scaffolds/contigs that are primarily microbial 
    (e.g., Wolbachia embedded in an insect genome). Generates a blacklist for downstream warnings.
    """
    input:
        genome=lambda wildcards: get_decontam_reference("host_genome", config["reference"]["genome"])
    output:
        kraken_out=temp(f"{DECONTAM_REF_DIR}/reference_audit.kraken.out"),
        report=f"{DECONTAM_REF_DIR}/reference_audit.kraken.report.tsv",
        blacklist=f"{DECONTAM_REF_DIR}/contam_scaffolds_blacklist.tsv"
    log:
        "logs/decontam/audit_reference_genome.log"
    benchmark:
        "benchmarks/decontam/audit_reference_genome.tsv"
    conda:
        "../../envs/decontam_kraken2.yaml"
    threads: CLASSIFIER_THREADS
    resources:
        mem_mb=8000
    params:
        db=CLASSIFIER_DB,
        nontarget_taxids=",".join(map(str, DECONTAM_CONFIG.get("policy", {}).get("nontarget_taxids", []))),
        uncertain_taxids=",".join(map(str, DECONTAM_CONFIG.get("policy", {}).get("uncertain_taxids", [])))
    shell:
        """
        mkdir -p {DECONTAM_REF_DIR} $(dirname {log}) benchmarks/decontam
        
        if [ -n "{params.db}" ] && [ -f "{params.db}/hash.k2d" ]; then
            k2 classify \
                --db {params.db} \
                --threads {threads} \
                --use-names \
                --output {output.kraken_out} \
                --report {output.report} \
                {input.genome} >> {log} 2>&1
                
            TAXIDS="{params.nontarget_taxids},{params.uncertain_taxids}"
            # Remove leading/trailing commas if any
            TAXIDS=$(echo "$TAXIDS" | sed 's/^,//;s/,$//;s/,,/,/g')
            
            echo -e "Scaffold_ID\tLength\tTaxID\tContam_Status" > {output.blacklist}
            if [ -n "$TAXIDS" ]; then
                taxonkit list --data-dir {params.db}/taxonomy --ids "$TAXIDS" > {output.kraken_out}.taxids 2>>{log} || touch {output.kraken_out}.taxids
                awk 'FNR==NR {{bad[$1]; next}} {{if ( ($3 in bad) ) print $2"\\t"$4"\\t"$3"\\tContam"}}' {output.kraken_out}.taxids {output.kraken_out} >> {output.blacklist}
                rm -f {output.kraken_out}.taxids
            fi
        else
            touch {output.kraken_out} {output.report} {output.blacklist}
            echo "Scaffold_ID\tLength\tTaxID\tContam_Status" > {output.blacklist}
        fi
        printf "[decontam_audit_reference_genome] audited reference genome\n" >> {log}
        """


rule decontam_prepare_host_reference:
    """
    Build the consolidated host rescue reference for Bowtie2.
    Transcriptome is concatenated before genome so the rescue layer stays transcript-aware.
    """
    input:
        transcriptome=lambda wildcards: get_decontam_reference("host_transcriptome", config["reference"]["transcriptome"]),
        genome=lambda wildcards: get_decontam_reference("host_genome", config["reference"]["genome"])
    output:
        reference=f"{DECONTAM_REF_DIR}/host_rescue.fa"
    log:
        "logs/decontam/host_reference.log"
    benchmark:
        "benchmarks/decontam/host_reference.tsv"
    shell:
        """
        mkdir -p {DECONTAM_REF_DIR} $(dirname {log}) benchmarks/decontam
        cat {input.transcriptome} {input.genome} > {output.reference}
        printf "[decontam_prepare_host_reference] built consolidated host rescue reference\n" > {log}
        """


rule decontam_build_host_bowtie2_index:
    """
    Build Bowtie2 index for host rescue automatically when absent.
    """
    input:
        reference=f"{DECONTAM_REF_DIR}/host_rescue.fa"
    output:
        index_dir=directory(f"{DECONTAM_INDEX_DIR}/bowtie2/host_rescue")
    log:
        "logs/decontam/host_bowtie2_index.log"
    benchmark:
        "benchmarks/decontam/host_bowtie2_index.tsv"
    conda:
        "../../envs/decontam_bowtie2.yaml"
    threads: 8
    resources:
        mem_mb=16000
    params:
        index_prefix=f"{DECONTAM_INDEX_DIR}/bowtie2/host_rescue/host_rescue"
    shell:
        """
        mkdir -p {output.index_dir} $(dirname {log}) benchmarks/decontam
        bowtie2-build --threads {threads} {input.reference} {params.index_prefix} > {log} 2>&1
        """


rule decontam_prepare_technical_reference:
    """
    Build the consolidated technical contamination reference for Bowtie2 prescreen.
    If no references are configured, emit a harmless dummy FASTA so the downstream collector can stay runnable.
    """
    input:
        references=lambda wildcards: get_reference_list("technical_contam")
    output:
        reference=f"{DECONTAM_REF_DIR}/technical_contam.fa"
    log:
        "logs/decontam/technical_reference.log"
    benchmark:
        "benchmarks/decontam/technical_reference.tsv"
    shell:
        """
        mkdir -p {DECONTAM_REF_DIR} $(dirname {log}) benchmarks/decontam
        if [ -n "{input.references}" ]; then
            cat {input.references} > {output.reference}
        else
            cat > {output.reference} <<'EOF'
>decontam_dummy
ACGTACGTACGTACGT
EOF
        fi
        printf "[decontam_prepare_technical_reference] built technical contamination reference\n" > {log}
        """


rule decontam_build_technical_bowtie2_index:
    """
    Build Bowtie2 index for technical contamination prescreen automatically when absent.
    """
    input:
        reference=f"{DECONTAM_REF_DIR}/technical_contam.fa"
    output:
        index_dir=directory(f"{DECONTAM_INDEX_DIR}/bowtie2/technical_contam")
    log:
        "logs/decontam/technical_bowtie2_index.log"
    benchmark:
        "benchmarks/decontam/technical_bowtie2_index.tsv"
    conda:
        "../../envs/decontam_bowtie2.yaml"
    threads: 8
    resources:
        mem_mb=16000
    params:
        index_prefix=f"{DECONTAM_INDEX_DIR}/bowtie2/technical_contam/technical_contam"
    shell:
        """
        mkdir -p {output.index_dir} $(dirname {log}) benchmarks/decontam
        bowtie2-build --threads {threads} {input.reference} {params.index_prefix} > {log} 2>&1
        """


rule decontam_prescreen:
    """
    Collect high-confidence technical contamination evidence.
    v1 keeps this as a placeholder ID collector so the pair-decision layer can already operate on a stable contract.
    """
    input:
        r1="results/02.trimmed_data/{sample}_R1_trimmed.fastq.gz",
        r2="results/02.trimmed_data/{sample}_R2_trimmed.fastq.gz",
        technical_index_dir=f"{DECONTAM_INDEX_DIR}/bowtie2/technical_contam",
        technical_reference=f"{DECONTAM_REF_DIR}/technical_contam.fa"
    output:
        tech_ids=temp(f"{DECONTAM_TMP_DIR}" + "/{sample}.technical_ids.txt.gz"),
        stats=f"{DECONTAM_STATS_DIR}" + "/{sample}_prescreen_stats.tsv"
    log:
        "logs/decontam/{sample}_prescreen.log"
    benchmark:
        "benchmarks/decontam/{sample}_prescreen.tsv"
    conda:
        "../../envs/decontam_bowtie2.yaml"
    threads: 8
    resources:
        mem_mb=16000
    params:
        technical_index_prefix=f"{DECONTAM_INDEX_DIR}/bowtie2/technical_contam/technical_contam",
        max_insert=PRESCREEN_MAX_INSERT,
        extra=PRESCREEN_EXTRA
    shell:
        """
        mkdir -p {DECONTAM_TMP_DIR} {DECONTAM_STATS_DIR} $(dirname {log}) benchmarks/decontam
        tmpdir=$(mktemp -d {DECONTAM_TMP_DIR}/{wildcards.sample}.prescreen.XXXXXX)
        trap 'rm -rf "$tmpdir"' EXIT

        bowtie2 \
            -x {params.technical_index_prefix} \
            -1 {input.r1} -2 {input.r2} \
            --al-conc-gz "$tmpdir/technical_mapped_%.fq.gz" \
            --no-mixed --no-discordant \
            --sensitive-local \
            -X {params.max_insert} \
            -p {threads} \
            {params.extra} \
            -S /dev/null >> {log} 2>&1

        if [ -f "$tmpdir/technical_mapped_1.fq.gz" ] && [ -f "$tmpdir/technical_mapped_2.fq.gz" ]; then
            seqkit seq -n "$tmpdir/technical_mapped_1.fq.gz" "$tmpdir/technical_mapped_2.fq.gz" | awk '{{print $1}}' | LC_ALL=C sort -u | gzip -c > {output.tech_ids}
        else
            echo "# read_id" | gzip -c > {output.tech_ids}
        fi

        matched_pairs=$(zcat {output.tech_ids} | awk '!/^#/' | wc -l)
        
        echo -e "sample\\tstage\\tmatched_pairs" > {output.stats}
        echo -e "{wildcards.sample}\\tprescreen\\t${{matched_pairs}}" >> {output.stats}
        printf "[decontam_prescreen] collected technical contamination evidence for %s\n" {wildcards.sample} > {log}
        """


rule decontam_host_rescue:
    """
    Collect host and ERCC rescue evidence using Bowtie2.
    This is the first real collector in v1: it writes host/ERCC ID sets and unresolved paired FASTQ for downstream classification.
    """
    input:
        r1="results/02.trimmed_data/{sample}_R1_trimmed.fastq.gz",
        r2="results/02.trimmed_data/{sample}_R2_trimmed.fastq.gz",
        technical_ids=f"{DECONTAM_TMP_DIR}" + "/{sample}.technical_ids.txt.gz",
        host_index_dir=f"{DECONTAM_INDEX_DIR}/bowtie2/host_rescue",
        host_reference=f"{DECONTAM_REF_DIR}/host_rescue.fa"
    output:
        host_ids=temp(f"{DECONTAM_TMP_DIR}" + "/{sample}.host_ids.txt.gz"),
        ercc_ids=temp(f"{DECONTAM_TMP_DIR}" + "/{sample}.ercc_ids.txt.gz"),
        unresolved_r1=temp(f"{DECONTAM_TMP_DIR}" + "/{sample}_R1_unresolved.fastq.gz"),
        unresolved_r2=temp(f"{DECONTAM_TMP_DIR}" + "/{sample}_R2_unresolved.fastq.gz"),
        stats=f"{DECONTAM_STATS_DIR}" + "/{sample}_host_rescue_stats.tsv"
    log:
        "logs/decontam/{sample}_host_rescue.log"
    benchmark:
        "benchmarks/decontam/{sample}_host_rescue.tsv"
    conda:
        "../../envs/decontam_bowtie2.yaml"
    threads: 8
    resources:
        mem_mb=32000
    params:
        host_index_prefix=f"{DECONTAM_INDEX_DIR}/bowtie2/host_rescue/host_rescue",
        ercc_reference=lambda wildcards: get_decontam_reference("ercc", ""),
        max_insert=HOST_RESCUE_MAX_INSERT,
        extra=HOST_RESCUE_EXTRA
    shell:
        """
        mkdir -p {DECONTAM_TMP_DIR} {DECONTAM_STATS_DIR} $(dirname {log}) benchmarks/decontam
        tmpdir=$(mktemp -d {DECONTAM_TMP_DIR}/{wildcards.sample}.host_rescue.XXXXXX)
        trap 'rm -rf "$tmpdir"' EXIT

        bowtie2 \
            -x {params.host_index_prefix} \
            -1 {input.r1} -2 {input.r2} \
            --no-mixed --no-discordant \
            --sensitive \
            -X {params.max_insert} \
            -p {threads} \
            {params.extra} 2>> {log} | \
        samtools view -e '[NM]<=5 && mapq>=10' -f 2 - | \
        cut -f1 | uniq -c | awk '$1==2 {{print $2}}' > "$tmpdir/host_ids.raw"
        
        gzip -c "$tmpdir/host_ids.raw" > {output.host_ids}
        
        seqkit grep -v -f "$tmpdir/host_ids.raw" {input.r1} -o {output.unresolved_r1}
        seqkit grep -v -f "$tmpdir/host_ids.raw" {input.r2} -o {output.unresolved_r2}

        if [ -n "{params.ercc_reference}" ] && [ -s "{params.ercc_reference}" ]; then
            ercc_prefix="$tmpdir/ercc_index"
            bowtie2-build --threads {threads} "{params.ercc_reference}" "$ercc_prefix" >> {log} 2>&1
            bowtie2 \
                -x "$ercc_prefix" \
                -1 {input.r1} -2 {input.r2} \
                --al-conc-gz "$tmpdir/ercc_mapped_%.fq.gz" \
                --no-mixed --no-discordant \
                --sensitive \
                -X {params.max_insert} \
                -p {threads} \
                -S /dev/null >> {log} 2>&1

            if [ -f "$tmpdir/ercc_mapped_1.fq.gz" ] && [ -f "$tmpdir/ercc_mapped_2.fq.gz" ]; then
                seqkit seq -n "$tmpdir/ercc_mapped_1.fq.gz" "$tmpdir/ercc_mapped_2.fq.gz" | awk '{{print $1}}' | LC_ALL=C sort -u | gzip -c > {output.ercc_ids}
            else
                echo "# read_id" | gzip -c > {output.ercc_ids}
            fi
        else
            echo "# read_id" | gzip -c > {output.ercc_ids}
        fi

        host_pairs=$(zcat {output.host_ids} | awk '!/^#/' | wc -l)
        ercc_pairs=$(zcat {output.ercc_ids} | awk '!/^#/' | wc -l)
        unresolved_pairs=$(echo $(($(zcat {output.unresolved_r1} | wc -l) / 4)))
        
        echo -e "sample\\tstage\\thost_pairs\\tercc_pairs\\tunresolved_pairs" > {output.stats}
        echo -e "{wildcards.sample}\\thost_rescue\\t${{host_pairs}}\\t${{ercc_pairs}}\\t${{unresolved_pairs}}" >> {output.stats}
        
        printf "[decontam_host_rescue] collected Bowtie2 rescue evidence for %s\n" {wildcards.sample} >> {log}
        """


rule decontam_classify_unresolved:
    """
    Collect classification-derived non-target and uncertain evidence.
    v1 intentionally keeps this collector as a valid empty-contract placeholder while the Bowtie2 rescue path is being debugged.
    """
    input:
        r1=f"{DECONTAM_TMP_DIR}" + "/{sample}_R1_unresolved.fastq.gz",
        r2=f"{DECONTAM_TMP_DIR}" + "/{sample}_R2_unresolved.fastq.gz",
        technical_ids=f"{DECONTAM_TMP_DIR}" + "/{sample}.technical_ids.txt.gz",
        host_ids=f"{DECONTAM_TMP_DIR}" + "/{sample}.host_ids.txt.gz",
        ercc_ids=f"{DECONTAM_TMP_DIR}" + "/{sample}.ercc_ids.txt.gz"
    output:
        kraken_output=f"{DECONTAM_TMP_DIR}" + "/{sample}.kraken.out",
        report=f"{DECONTAM_TMP_DIR}" + "/{sample}.kraken.report.tsv",
        nontarget_ids=f"{DECONTAM_TMP_DIR}" + "/{sample}.nontarget_ids.txt.gz",
        uncertain_ids=f"{DECONTAM_TMP_DIR}" + "/{sample}.uncertain_ids.txt.gz",
        stats=f"{DECONTAM_STATS_DIR}" + "/{sample}_classification_stats.tsv"
    log:
        "logs/decontam/{sample}_classification.log"
    benchmark:
        "benchmarks/decontam/{sample}_classification.tsv"
    conda:
        "../../envs/decontam_kraken2.yaml"
    threads: CLASSIFIER_THREADS
    resources:
        mem_mb=8000
    params:
        db=CLASSIFIER_DB,
        confidence=CLASSIFIER_CONFIG.get("confidence", 0.3),
        minimum_hit_groups=CLASSIFIER_CONFIG.get("minimum_hit_groups", 3),
        report_minimizer_data=CLASSIFIER_CONFIG.get("report_minimizer_data", True),
        memory_mapping=CLASSIFIER_CONFIG.get("memory_mapping", True),
        nontarget_taxids=repr(DECONTAM_CONFIG.get("policy", {}).get("nontarget_taxids", [])),
        uncertain_taxids=repr(DECONTAM_CONFIG.get("policy", {}).get("uncertain_taxids", [])),
        ecological_as_uncertain=DECONTAM_CONFIG.get("policy", {}).get("ecological_as_uncertain", True)
    shell:
        """
        mkdir -p {DECONTAM_TMP_DIR} {DECONTAM_STATS_DIR} $(dirname {log}) benchmarks/decontam
        if [ -n "{params.db}" ] && [ -f "{params.db}/hash.k2d" ] && [ -f "{params.db}/opts.k2d" ] && [ -f "{params.db}/taxo.k2d" ]; then
            extra_args=""
            if [ "{params.memory_mapping}" = "True" ]; then
                extra_args="$extra_args --memory-mapping"
            fi
            if [ "{params.report_minimizer_data}" = "True" ]; then
                extra_args="$extra_args --report-minimizer-data"
            fi

            k2 classify \
                --db {params.db} \
                --paired \
                --minimum-base-quality 20 \
                --use-names \
                --output {output.kraken_output} \
                --report {output.report} \
                --confidence {params.confidence} \
                --minimum-hit-groups {params.minimum_hit_groups} \
                --threads {threads} \
                $extra_args \
                {input.r1} {input.r2} >> {log} 2>&1
        else
            : > {output.kraken_output}
            : > {output.report}
        fi

        # Not using python script anymore. Use taxonkit + awk natively.
        
        NONTARGET=$(echo '{params.nontarget_taxids}' | sed 's/[][ '\'\"']//g')
        UNCERTAIN=$(echo '{params.uncertain_taxids}' | sed 's/[][ '\'\"']//g')
        
        if [ -n "$NONTARGET" ] && [ -n "{params.db}" ]; then
            taxonkit list --data-dir {params.db}/taxonomy --ids "$NONTARGET" > $tmpdir/nontarget.taxids 2>/dev/null || touch $tmpdir/nontarget.taxids
        else
            touch $tmpdir/nontarget.taxids
        fi
        
        if [ -n "$UNCERTAIN" ] && [ -n "{params.db}" ]; then
            taxonkit list --data-dir {params.db}/taxonomy --ids "$UNCERTAIN" > $tmpdir/uncertain.taxids 2>/dev/null || touch $tmpdir/uncertain.taxids
        else
            touch $tmpdir/uncertain.taxids
        fi

        # Process standard kraken output and categorize IDs directly
        awk '
            BEGIN {{FS="\\t"; OFS="\\t"}}
            FNR==NR {{ n[$1]; next }}  # nontarget taxids
            FNR==NR+1 {{ u[$1]; next }}  # uncertain taxids
            {{
                status=$1; qname=$2; taxid=$3
                
                # strip kraken suffix if present (/1 or /2)
                sub(/\\/[12]$/, "", qname)
                
                id_type = ""
                if (status == "C") {{
                    class_count++
                    if (taxid in n) {{
                        print qname > "'$tmpdir/nontarget_ids.txt'"
                        n_count++
                    }} else if (taxid in u && "{params.ecological_as_uncertain}" == "True") {{
                        print qname > "'$tmpdir/uncertain_ids.txt'"
                        u_count++
                    }}
                }} else {{
                    unclass_count++
                }}
            }}
            END {{
                print "sample\\tstage\\tclassified_pairs\\tunclassified_pairs\\tnontarget_pairs\\tuncertain_pairs" > "{output.stats}"
                printf "%s\\tclassification\\t%d\\t%d\\t%d\\t%d\\n", "{wildcards.sample}", class_count+0, unclass_count+0, n_count+0, u_count+0 > "{output.stats}"
            }}
        ' $tmpdir/nontarget.taxids $tmpdir/uncertain.taxids {output.kraken_output}
        
        [ -f $tmpdir/nontarget_ids.txt ] && sort -u $tmpdir/nontarget_ids.txt | gzip -c > {output.nontarget_ids} || echo "# read_id" | gzip -c > {output.nontarget_ids}
        [ -f $tmpdir/uncertain_ids.txt ] && sort -u $tmpdir/uncertain_ids.txt | gzip -c > {output.uncertain_ids} || echo "# read_id" | gzip -c > {output.uncertain_ids}

        printf "[decontam_classify_unresolved] collected k2 classification evidence for %s\n" {wildcards.sample} >> {log}
        """


rule decontam_pair_decision:
    """
    Single-pass pair-aware decision step.
    It loads evidence ID sets into memory and performs exactly one synchronized traversal over the trimmed FASTQ pair.
    """
    input:
        r1="results/02.trimmed_data/{sample}_R1_trimmed.fastq.gz",
        r2="results/02.trimmed_data/{sample}_R2_trimmed.fastq.gz",
        tech_ids=f"{DECONTAM_TMP_DIR}" + "/{sample}.technical_ids.txt.gz",
        host_ids=f"{DECONTAM_TMP_DIR}" + "/{sample}.host_ids.txt.gz",
        ercc_ids=f"{DECONTAM_TMP_DIR}" + "/{sample}.ercc_ids.txt.gz",
        nontarget_ids=f"{DECONTAM_TMP_DIR}" + "/{sample}.nontarget_ids.txt.gz",
        uncertain_ids=f"{DECONTAM_TMP_DIR}" + "/{sample}.uncertain_ids.txt.gz",
        prescreen_stats=f"{DECONTAM_STATS_DIR}" + "/{sample}_prescreen_stats.tsv",
        rescue_stats=f"{DECONTAM_STATS_DIR}" + "/{sample}_host_rescue_stats.tsv",
        classification_stats=f"{DECONTAM_STATS_DIR}" + "/{sample}_classification_stats.tsv"
    output:
        clean_r1=protected(f"{DECONTAM_CLEAN_DIR}" + "/{sample}_R1_clean.fastq.gz"),
        clean_r2=protected(f"{DECONTAM_CLEAN_DIR}" + "/{sample}_R2_clean.fastq.gz"),
        uncertain_r1=f"{DECONTAM_AUDIT_DIR}" + "/uncertain/{sample}_R1_uncertain.fastq.gz",
        uncertain_r2=f"{DECONTAM_AUDIT_DIR}" + "/uncertain/{sample}_R2_uncertain.fastq.gz",
        removed_r1=f"{DECONTAM_AUDIT_DIR}" + "/removed/{sample}_R1_removed.fastq.gz",
        removed_r2=f"{DECONTAM_AUDIT_DIR}" + "/removed/{sample}_R2_removed.fastq.gz",
        summary=f"{DECONTAM_STATS_DIR}" + "/{sample}_decision_summary.tsv"
    log:
        "logs/decontam/{sample}_pair_decision.log"
    benchmark:
        "benchmarks/decontam/{sample}_pair_decision.tsv"
    threads: 2
    resources:
        mem_mb=12000
    params:
        mode=DECONTAM_CONFIG.get("mode", "conservative"),
        retain_unclassified=DECONTAM_CONFIG.get("policy", {}).get("retain_unclassified", True),
        keep_audit_fastq=decontam_audit_enabled()
    conda:
        "../../envs/decontam.yaml"
    shell:
        """
        mkdir -p $(dirname {output.clean_r1}) $(dirname {output.uncertain_r1}) $(dirname {output.removed_r1}) $(dirname {output.summary})

        # 1. Determine read ID sets using awk
        # Priority: ERCC > Host > Tech > NonTarget > Uncertain
        # If a read ID string is in a higher priority, we ignore it in lower ones.
        
        # We need a clean fastq: ERCC + Host + (Uncertain if retain_unclassified=True)
        # We need to filter the FASTQ files using seqkit.
        
        # We use a single awk pass over all ID files to figure out the fate of EVERY read pair reported.
        echo -e "sample\\tRetained_Host\\tRetained_ERCC\\tRetained_Uncertain\\tRemoved_Tech\\tRemoved_NonTarget" > {output.summary}
        
        tmpdir=$(mktemp -d {DECONTAM_TMP_DIR}/{wildcards.sample}.pair_decision.XXXXXX)
        trap 'rm -rf "$tmpdir"' EXIT
        
        awk -v retain_unc="{params.retain_unclassified}" '
            /^#/ {{ next }}
            FNR==NR {{ fate[$1]="Retained_ERCC"; e++; next }}
            FNR==1 {{ file++ }}
            file==1 {{ if (!fate[$1]) {{ fate[$1]="Retained_Host"; h++ }}; next }}
            file==2 {{ if (!fate[$1]) {{ fate[$1]="Removed_Tech"; t++ }}; next }}
            file==3 {{ if (!fate[$1]) {{ fate[$1]="Removed_NonTarget"; n++ }}; next }}
            file==4 {{ if (!fate[$1]) {{ fate[$1]="Retained_Uncertain"; u++ }}; next }}
            END {{
                for (id in fate) {{
                    f = fate[id]
                    if (f == "Retained_ERCC" || f == "Retained_Host" || (f == "Retained_Uncertain" && retain_unc == "True")) {{
                        print id > "'$tmpdir'/clean_ids.txt"
                    }} else if (f == "Retained_Uncertain") {{
                        print id > "'$tmpdir'/uncertain_ids.txt"
                    }} else {{
                        print id > "'$tmpdir'/removed_ids.txt"
                    }}
                }}
                printf "%s\\t%d\\t%d\\t%d\\t%d\\t%d\\n", "{wildcards.sample}", h+0, e+0, u+0, t+0, n+0 > "{output.summary}"
            }}
        ' <(zcat -f {input.ercc_ids}) <(zcat -f {input.host_ids}) <(zcat -f {input.tech_ids}) <(zcat -f {input.nontarget_ids}) <(zcat -f {input.uncertain_ids})
        
        touch $tmpdir/clean_ids.txt $tmpdir/uncertain_ids.txt $tmpdir/removed_ids.txt
        
        # 2. Extract sequences using seqkit
        seqkit grep -f $tmpdir/clean_ids.txt {input.r1} -o {output.clean_r1}
        seqkit grep -f $tmpdir/clean_ids.txt {input.r2} -o {output.clean_r2}
        
        if [ "{params.keep_audit_fastq}" = "True" ]; then
            seqkit grep -f $tmpdir/uncertain_ids.txt {input.r1} -o {output.uncertain_r1}
            seqkit grep -f $tmpdir/uncertain_ids.txt {input.r2} -o {output.uncertain_r2}
            seqkit grep -f $tmpdir/removed_ids.txt {input.r1} -o {output.removed_r1}
            seqkit grep -f $tmpdir/removed_ids.txt {input.r2} -o {output.removed_r2}
        else
            echo "" | gzip -c > {output.uncertain_r1}
            echo "" | gzip -c > {output.uncertain_r2}
            echo "" | gzip -c > {output.removed_r1}
            echo "" | gzip -c > {output.removed_r2}
        fi
        
        printf "[decontam_pair_decision] completed pair-aware extraction for %s\n" {wildcards.sample} >> {log}
        """


rule clean_reads_fastqc:
    input:
        forward_reads=f"{DECONTAM_CLEAN_DIR}" + "/{sample}_R1_clean.fastq.gz",
        reverse_reads=f"{DECONTAM_CLEAN_DIR}" + "/{sample}_R2_clean.fastq.gz"
    output:
        forward_html=f"{DECONTAM_QC_DIR}" + "/{sample}_R1_clean_fastqc.html",
        reverse_html=f"{DECONTAM_QC_DIR}" + "/{sample}_R2_clean_fastqc.html",
        forward_zip=f"{DECONTAM_QC_DIR}" + "/{sample}_R1_clean_fastqc.zip",
        reverse_zip=f"{DECONTAM_QC_DIR}" + "/{sample}_R2_clean_fastqc.zip"
    threads:
        config["medium"]["t"]
    resources:
        mem_mb=config["fastqc"]["memory"]
    conda:
        "../../envs/qc.yaml"
    params:
        output_dir=DECONTAM_QC_DIR
    log:
        "logs/decontam/{sample}_clean_fastqc.log"
    shell:
        """
        mkdir -p {params.output_dir} $(dirname {log})
        fastqc -q -t {threads} --memory {resources.mem_mb} -o {params.output_dir} {input.forward_reads} {input.reverse_reads} >> {log} 2>&1
        """


rule decontam_project_summary:
    input:
        summaries=expand(f"{DECONTAM_STATS_DIR}" + "/{sample}_decision_summary.tsv", sample=SAMPLES)
    output:
        f"{DECONTAM_STATS_DIR}/project_decontam_summary.tsv"
    log:
        "logs/decontam/project_summary.log"
    benchmark:
        "benchmarks/decontam/project_summary.tsv"
    shell:
        """
        mkdir -p {DECONTAM_STATS_DIR} $(dirname {log}) benchmarks/decontam
        cat > {output} <<'EOF'
# id: 'decontam_fate_summary'
# section_name: 'Decontamination Fate'
# description: 'Read pair tracking through the conservative decontamination module.'
# plot_type: 'bargraph'
# pconfig:
#   id: 'decontam_bargraph'
#   title: 'Read Pairs Fate Allocation'
#   ylab: 'Number of Read Pairs'
Sample\tRetained_Host\tRetained_ERCC\tRetained_Uncertain\tRemoved_Tech\tRemoved_NonTarget
EOF
        for summary in {input.summaries}; do
            tail -n +2 "$summary" >> {output}
        done
        printf "[decontam_project_summary] aggregated %s samples\n" "$(echo {input.summaries} | wc -w)" > {log}
        """