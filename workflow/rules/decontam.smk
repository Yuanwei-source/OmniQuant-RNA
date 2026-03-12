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
            --very-sensitive-local \
            -X {params.max_insert} \
            -p {threads} \
            {params.extra} \
            -S /dev/null >> {log} 2>&1

        if [ -f "$tmpdir/technical_mapped_1.fq.gz" ] && [ -f "$tmpdir/technical_mapped_2.fq.gz" ]; then
            seqkit seq -n "$tmpdir/technical_mapped_1.fq.gz" "$tmpdir/technical_mapped_2.fq.gz" | awk '{{print $1}}' | LC_ALL=C sort -u | gzip -c > {output.tech_ids}
        else
            python - <<'PY'
import gzip
with gzip.open("{output.tech_ids}", "wt") as handle:
    handle.write("# read_id\\n")
PY
        fi

        python - <<'PY'
import gzip

with gzip.open("{output.tech_ids}", "rt") as handle:
    matched_pairs = sum(1 for line in handle if line.strip() and not line.startswith("#"))

with open("{output.stats}", "w") as handle:
    handle.write("sample\tstage\tmatched_pairs\\n")
    handle.write("{wildcards.sample}\tprescreen\t{{matched_pairs}}\\n".format(matched_pairs=matched_pairs))
PY
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
            --al-conc-gz "$tmpdir/host_mapped_%.fq.gz" \
            --un-conc-gz "$tmpdir/unresolved_%.fq.gz" \
            --no-mixed --no-discordant \
            --very-sensitive-local \
            -X {params.max_insert} \
            -p {threads} \
            {params.extra} \
            -S /dev/null >> {log} 2>&1

        if [ -f "$tmpdir/unresolved_1.fq.gz" ] && [ -f "$tmpdir/unresolved_2.fq.gz" ]; then
            mv "$tmpdir/unresolved_1.fq.gz" {output.unresolved_r1}
            mv "$tmpdir/unresolved_2.fq.gz" {output.unresolved_r2}
        else
            python - <<'PY'
import gzip
gzip.open("{output.unresolved_r1}", "wt").close()
gzip.open("{output.unresolved_r2}", "wt").close()
PY
        fi

        if [ -f "$tmpdir/host_mapped_1.fq.gz" ] && [ -f "$tmpdir/host_mapped_2.fq.gz" ]; then
            seqkit seq -n "$tmpdir/host_mapped_1.fq.gz" "$tmpdir/host_mapped_2.fq.gz" | awk '{{print $1}}' | LC_ALL=C sort -u | gzip -c > {output.host_ids}
        else
            python - <<'PY'
import gzip
with gzip.open("{output.host_ids}", "wt") as handle:
    handle.write("# read_id\\n")
PY
        fi

        if [ -n "{params.ercc_reference}" ] && [ -s "{params.ercc_reference}" ]; then
            ercc_prefix="$tmpdir/ercc_index"
            bowtie2-build --threads {threads} "{params.ercc_reference}" "$ercc_prefix" >> {log} 2>&1
            bowtie2 \
                -x "$ercc_prefix" \
                -1 {input.r1} -2 {input.r2} \
                --al-conc-gz "$tmpdir/ercc_mapped_%.fq.gz" \
                --no-mixed --no-discordant \
                --very-sensitive-local \
                -X {params.max_insert} \
                -p {threads} \
                -S /dev/null >> {log} 2>&1

            if [ -f "$tmpdir/ercc_mapped_1.fq.gz" ] && [ -f "$tmpdir/ercc_mapped_2.fq.gz" ]; then
                seqkit seq -n "$tmpdir/ercc_mapped_1.fq.gz" "$tmpdir/ercc_mapped_2.fq.gz" | awk '{{print $1}}' | LC_ALL=C sort -u | gzip -c > {output.ercc_ids}
            else
                python - <<'PY'
import gzip
with gzip.open("{output.ercc_ids}", "wt") as handle:
    handle.write("# read_id\\n")
PY
            fi
        else
            python - <<'PY'
import gzip
with gzip.open("{output.ercc_ids}", "wt") as handle:
    handle.write("# read_id\\n")
PY
        fi

        python - <<'PY'
import gzip

def count_ids(path):
    with gzip.open(path, "rt") as handle:
        return sum(1 for line in handle if line.strip() and not line.startswith("#"))

def count_fastq_records(path):
    with gzip.open(path, "rt") as handle:
        return sum(1 for _ in handle) // 4

host_pairs = count_ids("{output.host_ids}")
ercc_pairs = count_ids("{output.ercc_ids}")
unresolved_pairs = count_fastq_records("{output.unresolved_r1}")

with open("{output.stats}", "w") as handle:
    handle.write("sample\tstage\thost_pairs\tercc_pairs\tunresolved_pairs\\n")
    handle.write("{wildcards.sample}\thost_rescue\t{{host_pairs}}\t{{ercc_pairs}}\t{{unresolved_pairs}}\\n".format(
        host_pairs=host_pairs,
        ercc_pairs=ercc_pairs,
        unresolved_pairs=unresolved_pairs,
    ))
PY
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
        kraken_output=temp(f"{DECONTAM_TMP_DIR}" + "/{sample}.kraken.out"),
        report=temp(f"{DECONTAM_TMP_DIR}" + "/{sample}.kraken.report.tsv"),
        nontarget_ids=temp(f"{DECONTAM_TMP_DIR}" + "/{sample}.nontarget_ids.txt.gz"),
        uncertain_ids=temp(f"{DECONTAM_TMP_DIR}" + "/{sample}.uncertain_ids.txt.gz"),
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

        python workflow/scripts/parse_kraken_classification.py \
            --kraken-output {output.kraken_output} \
            --nodes-file {params.db}/taxo.k2d \
            --taxonomy-nodes {params.db}/taxonomy/nodes.dmp \
            --sample {wildcards.sample} \
            --nontarget-taxids '{params.nontarget_taxids}' \
            --uncertain-taxids '{params.uncertain_taxids}' \
            --ecological-as-uncertain {params.ecological_as_uncertain} \
            --output-nontarget {output.nontarget_ids} \
            --output-uncertain {output.uncertain_ids} \
            --output-stats {output.stats} >> {log} 2>&1

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
    script:
        "../scripts/decontam_pair_decision.py"


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