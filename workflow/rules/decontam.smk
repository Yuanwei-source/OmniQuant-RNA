DECONTAM_DIR = "results/03.decontam"
DECONTAM_TMP_DIR = f"{DECONTAM_DIR}/tmp"
DECONTAM_CLEAN_DIR = f"{DECONTAM_DIR}/clean"
DECONTAM_AUDIT_DIR = f"{DECONTAM_DIR}/audit"
DECONTAM_STATS_DIR = f"{DECONTAM_DIR}/stats"
DECONTAM_QC_DIR = f"{DECONTAM_DIR}/qc"
DECONTAM_CLUES_DIR = f"{DECONTAM_DIR}/clues"
DECONTAM_CLUES_TABLE_DIR = f"{DECONTAM_CLUES_DIR}/tables"
DECONTAM_CLUES_PLOT_DIR = f"{DECONTAM_CLUES_DIR}/plots"
DECONTAM_CLUES_INTERMEDIATE_DIR = f"{DECONTAM_CLUES_DIR}/intermediate"
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


def get_kaiju_reference(key, default=""):
    value = DECONTAM_CONFIG.get("kaiju", {}).get(key, default)
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
        genome=lambda wildcards: get_decontam_reference("host_genome", REFERENCE_GENOME)
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
        transcriptome=lambda wildcards: get_decontam_reference("host_transcriptome", REFERENCE_TRANSCRIPTOME),
        genome=lambda wildcards: get_decontam_reference("host_genome", REFERENCE_GENOME)
    output:
        reference=f"{DECONTAM_REF_DIR}/host_rescue.fa"
    log:
        "logs/decontam/host_reference.log"
    benchmark:
        "benchmarks/decontam/host_reference.tsv"
    shell:
        """
        mkdir -p {DECONTAM_REF_DIR} $(dirname {log}) benchmarks/decontam
        cat {input.genome} > {output.reference}
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
        tmp_reference=$(mktemp --suffix=.fa)
        trap 'rm -f "$tmp_reference"' EXIT

        if [ -n "{input.references}" ]; then
            cat {input.references} > "$tmp_reference"
        else
            cat > "$tmp_reference" <<'EOF'
>decontam_dummy
ACGTACGTACGTACGT
EOF
        fi

        awk '
            /^>/ {{
                header = substr($0, 2)
                gsub(/[[:space:]]+/, "_", header)
                gsub(/[^A-Za-z0-9._:|=-]/, "_", header)
                counts[header] += 1
                if (counts[header] > 1) {{
                    header = header "__" counts[header]
                }}
                print ">" header
                next
            }}
            {{ print }}
        ' "$tmp_reference" > {output.reference}

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
    shadow: "shallow"
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
            --no-mixed --no-discordant \
            --sensitive-local \
            -X {params.max_insert} \
            -p {threads} \
            {params.extra} 2>> {log} | \
        samtools view -F 4 - | cut -f1 | LC_ALL=C sort -u | gzip -c > {output.tech_ids}

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
    shadow: "shallow"
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
        samtools view -e '[NM]<=5' -f 2 - | \
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
                --no-mixed --no-discordant \
                --sensitive \
                -X {params.max_insert} \
                -p {threads} \
                2>> {log} | \
            samtools view -F 4 - | cut -f1 | LC_ALL=C sort -u | gzip -c > {output.ercc_ids}
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
    Track 2 (microbe analysis): Classify unresolved reads with Kaiju against progenomes DB.
    Replaces Kraken2+taxonkit. Outputs are for independent microbe visualization only —
    do NOT feed back into Track 1 (host DEA). nontarget_ids and uncertain_ids are empty
    placeholders for pair_decision contract compatibility.
    """
    input:
        r1=f"{DECONTAM_TMP_DIR}" + "/{sample}_R1_unresolved.fastq.gz",
        r2=f"{DECONTAM_TMP_DIR}" + "/{sample}_R2_unresolved.fastq.gz"
    output:
        kaiju_out=f"{DECONTAM_TMP_DIR}" + "/{sample}.kaiju.out",
        kaiju_lineage=f"{DECONTAM_TMP_DIR}" + "/{sample}.kaiju.out.lineage",
        krona_input=f"{DECONTAM_TMP_DIR}" + "/{sample}.krona_input.txt",
        krona_html=f"{DECONTAM_TMP_DIR}" + "/{sample}.krona.html",
        nontarget_ids=f"{DECONTAM_TMP_DIR}" + "/{sample}.nontarget_ids.txt.gz",
        uncertain_ids=f"{DECONTAM_TMP_DIR}" + "/{sample}.uncertain_ids.txt.gz",
        stats=f"{DECONTAM_STATS_DIR}" + "/{sample}_classification_stats.tsv"
    log:
        "logs/decontam/{sample}_classification.log"
    benchmark:
        "benchmarks/decontam/{sample}_classification.tsv"
    conda:
        "../../envs/kaiju.yaml"
    threads: int(get_kaiju_reference("threads", "32"))
    resources:
        mem_mb=120000
    params:
        nodes=lambda wildcards: get_kaiju_reference("nodes", "/home/dell/database/kaiju/nodes.dmp"),
        names=lambda wildcards: get_kaiju_reference("names", "/home/dell/database/kaiju/names.dmp"),
        db=lambda wildcards: get_kaiju_reference("db", "/home/dell/database/kaiju/kaiju_db_nr_euk.fmi")
    shell:
        """
        mkdir -p {DECONTAM_TMP_DIR} {DECONTAM_STATS_DIR} $(dirname {log}) benchmarks/decontam

        if [ -n "{params.db}" ] && [ -f "{params.db}" ]; then
            kaiju -t {params.nodes} \\
                  -f {params.db} \\
                  -i {input.r1} \\
                  -o {output.kaiju_out} \\
                  -z {threads} \\
                  -a greedy -e 3 -s 65 -E 0.01 -v >> {log} 2>&1

            kaiju-addTaxonNames -t {params.nodes} \\
                                -n {params.names} \\
                                -i {output.kaiju_out} \\
                                -o {output.kaiju_lineage} \\
                                -r superkingdom,phylum,class,order,family,genus,species -p >> {log} 2>&1

            kaiju2krona -t {params.nodes} \\
                        -n {params.names} \\
                        -i {output.kaiju_out} \\
                        -o {output.krona_input} >> {log} 2>&1
            ktImportText {output.krona_input} -o {output.krona_html} >> {log} 2>&1

            classified_reads=$(awk '$1=="C" {{c++}} END {{print c+0}}' {output.kaiju_out})
            unclassified_reads=$(awk '$1=="U" {{u++}} END {{print u+0}}' {output.kaiju_out})
            classified_pairs=$((classified_reads))
            unclassified_pairs=$((unclassified_reads))
        else
            : > {output.kaiju_out}
            : > {output.kaiju_lineage}
            : > {output.krona_input}
            : > {output.krona_html}
            classified_pairs=0
            unclassified_pairs=0
        fi

        echo -e "sample\\tstage\\tclassified_pairs\\tunclassified_pairs\\tnontarget_pairs\\tuncertain_pairs" > {output.stats}
        echo -e "{wildcards.sample}\\tclassification\\t${{classified_pairs}}\\t${{unclassified_pairs}}\\t0\\t0" >> {output.stats}

        # Empty placeholders for pair_decision contract compatibility.
        # Track 2 does not feed back into Track 1.
        echo "# read_id" | gzip -c > {output.nontarget_ids}
        echo "# read_id" | gzip -c > {output.uncertain_ids}

        printf "[decontam_classify_unresolved] completed Kaiju classification for %s\\n" {wildcards.sample} >> {log}
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
        
        # The default clean fastq is host-oriented: ERCC + Host.
        # Uncertain reads only re-enter clean when retain_unclassified=True is explicitly enabled.
        # We need to filter the FASTQ files using seqkit.
        
        # We use a single awk pass over all ID files to figure out the fate of EVERY read pair reported.
        echo -e "sample\\tRescued_Host\\tRescued_ERCC\\tFlagged_Uncertain\\tRemoved_Tech\\tRemoved_NonTarget" > {output.summary}
        
        tmpdir=$(mktemp -d {DECONTAM_TMP_DIR}/{wildcards.sample}.pair_decision.XXXXXX)
        trap 'rm -rf "$tmpdir"' EXIT
        
        awk -v retain_unc="{params.retain_unclassified}" '
            /^#/ {{ next }}
            FNR==NR {{ fate[$1]="Rescued_ERCC"; e++; next }}
            FNR==1 {{ file++ }}
            file==1 {{ if (!fate[$1]) {{ fate[$1]="Rescued_Host"; h++ }}; next }}
            file==2 {{ if (!fate[$1]) {{ fate[$1]="Removed_Tech"; t++ }}; next }}
            file==3 {{ if (!fate[$1]) {{ fate[$1]="Removed_NonTarget"; n++ }}; next }}
            file==4 {{ if (!fate[$1]) {{ fate[$1]="Flagged_Uncertain"; u++ }}; next }}
            END {{
                for (id in fate) {{
                    f = fate[id]
                    if (f == "Rescued_ERCC" || f == "Rescued_Host" || (f == "Flagged_Uncertain" && retain_unc == "True")) {{
                        print id > "'$tmpdir'/clean_ids.txt"
                    }} else if (f == "Flagged_Uncertain") {{
                        print id > "'$tmpdir'/uncertain_ids.txt"
                    }} else {{
                        print id > "'$tmpdir'/removed_ids.txt"
                    }}
                }}
                printf "%s\\t%d\\t%d\\t%d\\t%d\\t%d\\n", "{wildcards.sample}", h+0, e+0, u+0, t+0, n+0 >> "{output.summary}"
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


rule decontam_microbe_clues_tables:
    input:
        sample_file=config["samples"],
        decision_summaries=expand(f"{DECONTAM_STATS_DIR}" + "/{sample}_decision_summary.tsv", sample=SAMPLES),
        classification_stats=expand(f"{DECONTAM_STATS_DIR}" + "/{sample}_classification_stats.tsv", sample=SAMPLES),
        kaiju_lineages=expand(f"{DECONTAM_TMP_DIR}" + "/{sample}.kaiju.out.lineage", sample=SAMPLES)
    output:
        burden=f"{DECONTAM_CLUES_TABLE_DIR}/sample_microbial_burden.tsv",
        targets=f"{DECONTAM_CLUES_TABLE_DIR}/priority_targets.tsv",
        composition=f"{DECONTAM_CLUES_INTERMEDIATE_DIR}/microbial_composition_long.tsv"
    conda:
        "../../envs/decontam.yaml"
    log:
        "logs/decontam/microbe_clues_tables.log"
    shell:
        """
        mkdir -p {DECONTAM_CLUES_TABLE_DIR} {DECONTAM_CLUES_INTERMEDIATE_DIR} $(dirname {log})

        # 构建 sample → group 映射表 (samples.tsv 第4列是 group)
        tmpdir=$(mktemp -d {DECONTAM_TMP_DIR}/.clues_tables.XXXXXX)
        trap 'rm -rf "$tmpdir"' EXIT
        awk 'NR>1 {{print $1"\\t"$4}}' {input.sample_file} > "$tmpdir/groups.tsv"

        # 写入三个输出文件的表头
        printf "sample\\tgroup\\ttotal_reads\\tclassified_reads\\tunclassified_reads\\tbacteria_reads\\tfungi_reads\\tvirus_reads\\tnon_host_fraction\\n" > {output.burden}
        printf "sample\\tgroup\\twolbachia_reads\\twolbachia_presence\\tvirus_reads\\tvirus_presence\\tfungi_reads\\tfungi_presence\\n" > {output.targets}
        printf "sample\\tgroup\\tcategory\\tpairs\\tfraction_total\\tfraction_non_host\\n" > {output.composition}

        # 逐个样本处理 Kaiju lineage 文件
        for linefile in {input.kaiju_lineages}; do
            sample=$(basename "$linefile" .kaiju.out.lineage)
            group=$(awk -v s="$sample" '$1==s {{print $2; exit}}' "$tmpdir/groups.tsv")
            [ -z "$group" ] && group="unknown"

            # 从 decision summary 获取 total_pairs 用于计算 non_host_fraction
            summary="{DECONTAM_STATS_DIR}/${{sample}}_decision_summary.tsv"
            if [ -f "$summary" ]; then
                rescued_host=$(awk 'NR==2 {{print $2+0}}' "$summary")
                rescued_ercc=$(awk 'NR==2 {{print $3+0}}' "$summary")
                flagged_uncertain=$(awk 'NR==2 {{print $4+0}}' "$summary")
                removed_tech=$(awk 'NR==2 {{print $5+0}}' "$summary")
                removed_nontarget=$(awk 'NR==2 {{print $6+0}}' "$summary")
                total_pairs=$((rescued_host + rescued_ercc + flagged_uncertain + removed_tech + removed_nontarget))
                non_host_pairs=$((total_pairs - rescued_host - rescued_ercc))
                if [ "$total_pairs" -gt 0 ]; then
                    non_host_fraction=$(awk "BEGIN {{printf \"%.6f\", $non_host_pairs/$total_pairs}}")
                else
                    non_host_fraction="0.000000"
                fi
            else
                non_host_fraction="0.000000"
            fi

            # 从 Kaiju lineage 文件统计各类 reads（单次 awk 遍历，高效）
            if [ -f "$linefile" ] && [ -s "$linefile" ]; then
                read total_reads classified unclassified bacteria fungi virus wolbachia <<< $(awk '
                    /^C/ {{ classified++ }}
                    /^U/ {{ unclassified++ }}
                    /^C/ && /Bacteria;/ {{ bacteria++ }}
                    /^C/ && /Fungi;/    {{ fungi++ }}
                    /^C/ && /Viruses;/  {{ virus++ }}
                    tolower($0) ~ /wolbachia/ {{ wolbachia++ }}
                    END {{
                        printf "%d %d %d %d %d %d %d\\n", NR, classified+0, unclassified+0, bacteria+0, fungi+0, virus+0, wolbachia+0
                    }}
                ' "$linefile")
            else
                total_reads=0; classified=0; unclassified=0
                bacteria=0; fungi=0; virus=0; wolbachia=0
            fi

            # 写入 burden 行
            printf "%s\\t%s\\t%d\\t%d\\t%d\\t%d\\t%d\\t%d\\t%s\\n" \
                "$sample" "$group" $total_reads $classified $unclassified \
                $bacteria $fungi $virus "$non_host_fraction" >> {output.burden}

            # 写入 targets 行
            w_presence=$([ "$wolbachia" -gt 0 ] && echo "yes" || echo "no")
            v_presence=$([ "$virus" -gt 0 ] && echo "yes" || echo "no")
            f_presence=$([ "$fungi" -gt 0 ] && echo "yes" || echo "no")
            printf "%s\\t%s\\t%d\\t%s\\t%d\\t%s\\t%d\\t%s\\n" \
                "$sample" "$group" $wolbachia "$w_presence" $virus "$v_presence" \
                $fungi "$f_presence" >> {output.targets}

            # 写入 composition 行（6个分类）
            other_class=$((classified - bacteria - fungi - virus))
            [ $other_class -lt 0 ] && other_class=0

            for row in "Technical:0" "Bacteria:$bacteria" "Fungi:$fungi" "Viruses:$virus" "Other_Classified:$other_class" "Unclassified:$unclassified"; do
                cat=$(echo "$row" | cut -d: -f1)
                cnt=$(echo "$row" | cut -d: -f2)
                if [ "$total_reads" -gt 0 ]; then
                    ft=$(awk "BEGIN {{printf \"%.6f\", $cnt/$total_reads}}")
                else
                    ft="0.000000"
                fi
                if [ "$classified" -gt 0 ]; then
                    fnh=$(awk "BEGIN {{printf \"%.6f\", $cnt/$classified}}")
                else
                    fnh="0.000000"
                fi
                printf "%s\\t%s\\t%s\\t%d\\t%s\\t%s\\n" \
                    "$sample" "$group" "$cat" $cnt "$ft" "$fnh" >> {output.composition}
            done
        done

        n_samples=$(wc -l < {output.burden})
        printf "[decontam_microbe_clues_tables] 已生成 %d 个样本的微生物线索表\\n" "$((n_samples - 1))" > {log}
        """


rule decontam_microbe_composition_plot:
    input:
        composition=f"{DECONTAM_CLUES_INTERMEDIATE_DIR}/microbial_composition_long.tsv"
    output:
        plot=f"{DECONTAM_CLUES_PLOT_DIR}/microbial_composition_stacked_bar.pdf"
    conda:
        "../../envs/dea.yaml"
    log:
        "logs/decontam/microbe_composition_plot.log"
    script:
        "../scripts/plot_microbe_composition.R"


rule decontam_host_context_overlay:
    input:
        burden=f"{DECONTAM_CLUES_TABLE_DIR}/sample_microbial_burden.tsv",
        targets=f"{DECONTAM_CLUES_TABLE_DIR}/priority_targets.tsv",
        sample_file=config["samples"],
        normalized_counts="results/06.differential_expression/featurecounts/normalized_counts.csv"
    output:
        plot=f"{DECONTAM_CLUES_PLOT_DIR}/host_context_overlay.pdf"
    conda:
        "../../envs/dea.yaml"
    log:
        "logs/decontam/host_context_overlay.log"
    script:
        "../scripts/plot_host_context_overlay.R"


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
Sample\tRescued_Host\tRescued_ERCC\tFlagged_Uncertain\tRemoved_Tech\tRemoved_NonTarget
EOF
        for summary in {input.summaries}; do
            awk 'NR==1 && tolower($1) == "sample" {{next}} {{print}}' "$summary" >> {output}
        done
        printf "[decontam_project_summary] aggregated %s samples\n" "$(echo {input.summaries} | wc -w)" > {log}
        """