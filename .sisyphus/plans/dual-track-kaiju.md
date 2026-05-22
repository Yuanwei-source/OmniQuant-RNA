# 双轨流水线重构计划

## TL;DR

> **目标**: 将去污染模块拆分为独立的 Track 1（宿主纯净表达流）和 Track 2（微生物互作流），Track 2 用 Kaiju 氨基酸比对替代 Kraken2 + taxonkit。
>
> **核心收益**: Track 2 从 132GB DB + 240 万 taxid 展开降为 Kaiju nr_euk（~40GB）+ grep lineage，同时提升对高突变共生菌/RNA 病毒的敏感度。Track 1 不受影响。
>
> **预计改动**: 1 个新 conda env + 1 条规则重写 + 2 个新脚本 + 线索输出格式更新。

---

## Context

### 当前问题

| 问题 | 详情 |
|------|------|
| Kraken2 太重 | 132GB 全量 nt DB，对低污染样本浪费资源 |
| taxonkit 脆弱 | 展开 240 万 taxid，读格式依赖缩进，跨 DB 版本易失效 |
| DNA 层面敏感度低 | 对高 AT/突变共生菌、RNA 病毒不友好 |
| 匹配逻辑复杂 | awk 多文件关联 + taxonkit 静默失败难以调试 |
| 对共识 DEA 零影响 | 三次不同配置跑出来 Tier A 全是 1287 |

### 新设计原则

```
Track 1 (纯净表达流)              Track 2 (微生物互作流)
════════════════════════           ════════════════════════════
Raw FASTQ                          Track 1 的 unmapped reads
  │                                      │
fastp                                    ▼
  │                               Kaiju 翻译 6 框 → nr_euk
Bowtie2 宿主比对                         │
  ├─ mapped → Host Clean FASTQ     分类结果 (.out)
  │   ↓                                  │
  │   HISAT2                        kaiju-addTaxonNames
  │   4 quantifiers                       │
  │   DESeq2                        lineage 标注完成
  │   Consensus                           │
  │                               ┌───────┴───────┐
  └─ unmapped ──────────────────→  │ grep lineage   │
                                   │ 靶向提取       │
                                   └───────┬───────┘
                                           │
                                    seqtk subseq → Symbiont.fq
                                           │
                                    ┌──────┴──────┐
                                    │ Krona HTML  │  MEGAHIT/Trinity
                                    │ 交互饼图    │  从头组装(可选)
                                    └─────────────┘
```

**Track 2 不回灌 Track 1**——维持现有设计原则。

---

## Work Objectives

### Core Objective
用 Kaiju 氨基酸分类管线替代 Kraken2 + taxonkit，实现轻量、精准、可审计的微生物线索侧支。

### Concrete Deliverables

| 文件 | 说明 |
|------|------|
| `envs/kaiju.yaml` | Kaiju + kaiju-addTaxonNames + Krona conda 环境 |
| `workflow/scripts/extract_symbiont_reads.sh` | grep lineage + seqtk 靶向提取 |
| `workflow/scripts/build_krona_report.sh` | kaiju2krona → Krona HTML |
| `workflow/rules/decontam.smk` | 重写 `decontam_classify_unresolved` 及 clues 相关规则 |
| `config/config.yaml` | 新增 `decontam.kaiju` 配置节 |
| `results/03.decontam/clues/` | 更新输出格式 |

---

## TODOs

- [x] 1. **搭建 Kaiju conda 环境**

  **What to do**:
  - 创建 `envs/kaiju.yaml`：
    ```yaml
    name: kaiju
    channels:
      - conda-forge
      - bioconda
    dependencies:
      - kaiju>=1.9
      - krona>=2.8
      - seqtk>=1.4
      - seqkit>=2.5
    ```
  - `conda env create -f envs/kaiju.yaml`
  - 确认 Kaiju 可运行: `kaiju -h`
  - 准备数据库：下载 nr_euk（~40GB）或使用已有的 nr 蛋白库
    ```bash
    # Kaiju 官方创建 DB 命令
    kaiju-makedb -s nr_euk -t 32
    ```
  - DB 路径放入 config: `/mnt/nas/Database/kaiju_nr_euk`

  **QA**:
  ```
  Scenario: Kaiju 能对测试 reads 分类
    Tool: bash
    Steps: echo ">test\nACGT" | kaiju -t nodes.dmp -f nr_euk.fmi -i /dev/stdin
    Expected: 正常退出，无报错
  ```

  **Commit**: NO

- [x] 2. **重写 decontam_classify_unresolved 规则**

  **What to do**:
  - 修改 `workflow/rules/decontam.smk` 中的 `decontam_classify_unresolved` 规则
  - **移除**：Kraken2 + taxonkit 分类 + nontarget/uncertain taxid 匹配
  - **替换为**：Kaiju 三步管线
    ```bash
    # Step A: Kaiju 6-框翻译 + 蛋白比对
    kaiju -t {params.nodes} -f {params.db} \
      -i {input.unresolved_r1} -j {input.unresolved_r2} \
      -z {threads} -o {output.kaiju_out}

    # Step B: 加 lineage
    kaiju-addTaxonNames -t {params.nodes} -n {params.names} \
      -i {output.kaiju_out} -o {output.kaiju_out}.lineage \
      -r superkingdom,phylum,class,order,family,genus,species

    # Step C: Krona 图
    kaiju2krona -t {params.nodes} -n {params.names} \
      -i {output.kaiju_out} -o {output.krona_in}
    ktImportText {output.krona_in} -o {output.krona_html}
    ```
  - 输出文件：
    - `{sample}.kaiju.out` — 原始分类结果
    - `{sample}.kaiju.out.lineage` — 带 lineage 的标注版
    - `{sample}.krona.html` — Krona 交互饼图
  - 移除原有 `nontarget_ids`、`uncertain_ids`、`classification_stats` 输出
  - 这些原先传给 `pair_decision` 的文件 → 改为空占位（Track 2 不回灌 Track 1）

  **Must NOT do**:
  - 不要改 `pair_decision` 规则——它对 nontarget/uncertain 空输入已有容错
  - 不要改宿主救援规则

  **QA**:
  ```
  Scenario: Kaiju 规则产生 lineage 文件
    Tool: bash
    Steps: head -5 results/03.decontam/tmp/ND1.kaiju.out.lineage
    Expected: 每行含 domain;phylum;class;... 格式的 lineage
    Evidence: .sisyphus/evidence/kaiju-output-check.txt
  ```

  **Commit**: YES — `refactor: replace kraken2+taxonkit with kaiju in decontam classify rule`

- [x] 3. **创建靶向提取脚本**

  **What to do**:
  - 创建 `workflow/scripts/extract_symbiont_reads.sh`
  - 功能：从 Kaiju lineage 输出中 grep 指定分类群，用 seqtk 从 unmapped FASTQ 提取
  - 用法：
    ```bash
    extract_symbiont_reads.sh \
      --kaiju-lineage results/03.decontam/tmp/ND1.kaiju.out.lineage \
      --r1 results/03.decontam/tmp/ND1_R1_unresolved.fastq.gz \
      --r2 results/03.decontam/tmp/ND1_R2_unresolved.fastq.gz \
      --targets "Bacteria;Fungi;Viruses;" \
      --output-prefix results/03.decontam/clues/extracted/ND1
    ```
  - 逻辑：
    ```bash
    # 1. grep lineage 获取目标 reads ID
    grep -E "$TARGETS" "$LINEAGE" | cut -f2 > "$TMP/ids.txt"
    # 2. seqtk 提取
    seqtk subseq "$R1" "$TMP/ids.txt" | gzip > "${OUT}_symbiont_R1.fq.gz"
    seqtk subseq "$R2" "$TMP/ids.txt" | gzip > "${OUT}_symbiont_R2.fq.gz"
    ```
  - 默认 targets: `"Bacteria;"` `"Fungi;"` `"Viruses;"` `"Wolbachia;"`

  **QA**:
  ```
  Scenario: 能从芫菁数据提取细菌 reads
    Tool: bash
    Steps: ./extract_symbiont_reads.sh ... --targets "Xanthomonas;Escherichia;"
    Expected: 生成 *_symbiont_R1.fq.gz 且非空
    Evidence: .sisyphus/evidence/extract-symbiont-check.txt
  ```

  **Commit**: YES — `feat: add symbiont read extraction script using kaiju lineage + seqtk`

- [x] 4. **更新 clues 侧支汇总规则**

  **What to do**:
  - 当前的 `decontam_microbe_clues_tables` 和 `decontam_project_summary` 规则依赖 Kraken2 的 `report.tsv` 格式
  - **改为**读取 Kaiju 的 `.kaiju.out` 格式
  - Kaiju 输出格式（TSV）：
    ```
    C/V    read_id    taxid    score    ...
    ```
  - `summarize_microbe_clues.py` 需改为解析此格式，按 lineage 层级聚合
  - 替代方案：直接从 Kaiju lineage 输出统计各分类层级 reads 数（更简单）
    ```bash
    grep -c "Bacteria;" results/03.decontam/tmp/ND1.kaiju.out.lineage
    grep -c "Fungi;" results/03.decontam/tmp/ND1.kaiju.out.lineage
    grep -c "Viruses;" results/03.decontam/tmp/ND1.kaiju.out.lineage
    ```
  - `priority_targets.tsv`：grep 特定属名（如 "Wolbachia;"）
  
  **QA**:
  ```
  Scenario: clues/sample_microbial_burden.tsv 有非零细菌比例
    Tool: bash
    Steps: grep -v "^sample" results/03.decontam/clues/tables/sample_microbial_burden.tsv | awk -F'\t' '{print $1,$13}'
    Expected: 每样本的 bacteria_pairs > 0
    Evidence: .sisyphus/evidence/clues-updated-check.txt
  ```

  **Commit**: YES — `refactor: update microbe clues tables for kaiju lineage format`

- [x] 5. **新增 config 配置节**

  **What to do**:
  - 在 `config/config.yaml` 的 `decontam` 下新增 `kaiju` 节：
    ```yaml
    decontam:
      kaiju:
        enabled: true
        db: "/mnt/nas/Database/kaiju_nr_euk/nr_euk.fmi"
        nodes: "/mnt/nas/Database/kaiju_nr_euk/nodes.dmp"
        names: "/mnt/nas/Database/kaiju_nr_euk/names.dmp"
        threads: 32
      extract:
        targets: ["Bacteria;", "Fungi;", "Viruses;", "Wolbachia;"]
      clues:
        top_taxa_n: 3
        priority_targets:
          wolbachia: "Wolbachia;"
          virus: "Viruses;"
          fungi: "Fungi;"
          bacteria: "Bacteria;"
    ```
  - 移除/注释旧的 `classifier` 和 `policy` 节（Kraken2 相关）

  **QA**:
  ```
  Scenario: config 解析无误
    Tool: bash
    Steps: python3 -c "import yaml; c=yaml.safe_load(open('config/config.yaml')); print(c['decontam']['kaiju']['db'])"
    Expected: 打印 DB 路径
  ```

  **Commit**: YES — `feat: add kaiju config section, deprecate kraken2 classifier config`

- [ ] 6. **端到端测试（芫菁数据）**

  **What to do**:
  - 用新规则跑 ON 管道（Track 2 走 Kaiju）
  - 验证：
    - Track 1 输出不变（共识 Tier A = 1,287）
    - Kaiju lineage 输出有细菌/真菌/病毒分类
    - Krona HTML 可打开
    - `extract_symbiont_reads.sh` 能提取到目标 reads
    - clues 表更新正确

  **QA**:
  ```
  Scenario: Track 1 共识结果不变
    Tool: bash
    Steps: grep "tier_a_n" results.backup.epicauta_kaiju/07.consensus_expression/diapause_vs_non-diapause/consensus_summary.tsv
    Expected: tier_a_n = 1287

  Scenario: Track 2 有细菌分类
    Tool: bash
    Steps: grep -c "Bacteria;" results.backup.epicauta_kaiju/03.decontam/tmp/ND1.kaiju.out.lineage
    Expected: > 100
    Evidence: .sisyphus/evidence/kaiju-e2e-check.txt
  ```

  **Commit**: NO (results only)

- [x] 7. **更新论文 §2.4 和 §3.6**

  **What to do**:
  - Methods §2.4 (Decontamination Module)：
    - 描述双轨架构：Track 1 Bowtie2 宿主救援 → clean FASTQ；Track 2 Kaiju 氨基酸分类
    - 强调 Track 2 不回灌 Track 1
  - Results §3.6 (Epicauta)：
    - 报告 Kaiju 分类结果（细菌/真菌/病毒 reads 数）
    - 对比 PlusPFP / 全量 nt / Kaiju 三种分类策略
    - 论述 Kaiju 对非模式昆虫的优势（氨基酸比对 → 密码子不敏感）

  **Must NOT do**:
  - 不要说 "ground truth"
  - 不要声称 Kaiju 比 Kraken2 "更好"（说 "更适合非模式昆虫场景"）

  **Commit**: YES — `docs: update methods/results for dual-track kaiju pipeline`

---

## Execution Strategy

```
Task 1 (Kaiju env) ──┐
                      ├── Task 2 (重写规则) ── Task 3 (提取脚本) ── Task 4 (clues)
Task 5 (config) ──────┘                                              │
                                                          Task 6 (E2E测试)
                                                                      │
                                                          Task 7 (论文更新)
```

Task 1 和 5 可并行。Task 2-4 依赖 1。Task 6 依赖全部。Task 7 依赖 6。

---

## Verification Strategy

- Track 1 不受影响：共识 Tier A 结果与旧版一致（1,287）
- Track 2 产生 Kaiju lineage 文件，含非零细菌/真菌计数
- Krona HTML 可在浏览器中打开
- 提取脚本能产出 Symbiont.fq
- 论文更新后 §2.4 描述双轨、§3.6 有 Kaiju 对比数据

---

## Success Criteria

- [ ] Kaiju conda 环境创建成功，kaiju -h 正常
- [ ] 重跑 ON 管道后 consensus Tier A = 1,287（Track 1 不变）
- [ ] `*.kaiju.out.lineage` 含 Bacteria/Fungi/Viruses 分类
- [ ] `clues/sample_microbial_burden.tsv` 含非零 bacteria_pairs
- [ ] Krona HTML 可生成
- [ ] 论文 §2.4 和 §3.6 已更新
