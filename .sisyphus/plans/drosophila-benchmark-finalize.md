# Drosophila Benchmark 收尾计划

## TL;DR

> **目标**: 将分散在两个备份目录中的 benchmark 输出合并到 `results/benchmark/`，更新 `benchmark_summary.md` 纳入 decontam COMPARE 结果，执行最终一致性检查。
>
> **无需重跑管道**: ON/OFF 备份中已有完整的 pipeline 输出（stage 00-08 + benchmark）。
>
> **预计耗时**: 10 分钟（纯文件操作 + 文档更新 + grep 检查）。

---

## Context

### 当前状态

- `results.backup.decontam_on/benchmark/` — 有 Task 1/2/5/6/7 的所有输出（15 个文件）
- `results.backup.decontam_off/benchmark/` — 同样有（与 ON 相同，因为 ablation 依赖 ON 结果）
- `results/benchmark/` — 仅有 decontam COMPARE 输出（刚生成）
- `results/` — pipeline 仅跑到 stage 05，无 DEA/consensus
- `benchmark_decontam_spikein.R` — Mode 2/3 刚补全，COMPARE 已验证通过

### 为何不需要重跑

| 考虑 | 结论 |
|------|------|
| 备份结果是真实管道产出 | ✅ 来自 decontam ON/OFF 两遍完整运行 |
| 当前 `results/` 部分覆盖 | ⚠️ 仅 stage 00-05，但 benchmark 需要 stage 06-07 |
| 重跑管道耗时 | ❌ 数小时，且结果应与备份一致 |
| 直接从备份合并 | ✅ 结果等价，零时间成本 |

### Plan 文档标记问题

Plan 中所有 task 标记为 `[x]` 是因为**脚本已写完**，不是**结果已合并**。本次收尾将结果统一并做最终检查。

---

## Work Objectives

### Core Objective
将果蝇 benchmark 的所有输出合并到一个 `results/benchmark/` 目录下，更新摘要文档，完成一致性检查。

### Concrete Deliverables
- `results/benchmark/` 包含全部 6 类 benchmark 输出
- `results/benchmark/benchmark_summary.md` 更新版（含 decontam COMPARE 数据）
- Task 8 一致性检查通过

---

## Verification Strategy

- **文件完整性**: `results/benchmark/` 下所有预期文件存在
- **摘要完整性**: benchmark_summary.md 包含 5 段落（degradation / subsampling / decontam / ablation / limitations）
- **语言纪律**: grep 不含 "ground truth" / "gold standard"
- **decontam 数据**: Tier A retention 96.77%、logFC r=0.9932 等关键数字出现在摘要中

---

## TODOs

- [ ] 1. 合并 Benchmark 输出文件到 `results/benchmark/`

  **What to do**:
  - 从 `results.backup.decontam_on/benchmark/` 复制 Task 1/2/5/6 输出到 `results/benchmark/`
  - 保留已存在的 decontam COMPARE 文件（不覆盖）
  - 合并后的 `results/benchmark/` 应包含：

  | 文件 | 来源 | Task |
  |------|------|------|
  | `annotation_degradation.tsv` | ON backup | Task 1 |
  | `annotation_degradation_summary.tsv` | ON backup | Task 1 |
  | `subsampling_stability.tsv` | ON backup | Task 2 |
  | `subsampling_summary.tsv` | ON backup | Task 2 |
  | `ablation_summary.tsv` | ON backup | Task 5 |
  | `statistical_tests.tsv` | ON backup | Task 6 |
  | `bombyx_case_study.md` | ON backup | Task 0 |
  | `decontam_assessment.txt` | ON backup | Task 3 |
  | `figures/*.pdf` (8 张) | ON backup | Task 1/2 |
  | `decontam_comparison_summary.tsv` | 已有 | Task 3 |
  | `decontam_gene_level_comparison.tsv` | 已有 | Task 3 |
  | `tier_transition_matrix.tsv` | 已有 | Task 3 |
  | `decontam_efficiency.tsv` | 已有 | Task 3 |
  | `figures/logFC_ON_vs_OFF.pdf` | 已有 | Task 3 |
  | `figures/decontam_read_fate.pdf` | 已有 | Task 3 |
  | `figures/tier_transition_heatmap.pdf` | 已有 | Task 3 |

  **Must NOT do**:
  - 不要覆盖刚生成的 decontam COMPARE 文件
  - 不要直接 mv 备份目录（保留备份）

  **QA Scenarios**:
  ```
  Scenario: 所有预期文件存在
    Tool: bash
    Steps: ls results/benchmark/*.tsv results/benchmark/*.md results/benchmark/figures/*.pdf
    Expected: 至少 12 个 .tsv/.txt/.md + 11 张图
    Evidence: .sisyphus/evidence/benchmark-files-check.txt
  ```

  **Commit**: NO (output artifacts)

- [ ] 2. 更新 `benchmark_summary.md` 纳入 decontam COMPARE 数据

  **What to do**:
  - 读现有 `benchmark_summary.md`（从 ON backup 复制过来后）
  - 在 "Decontam" 段落补充以下 COMPARE 数据：

  | 指标 | 值 | 含义 |
  |------|-----|------|
  | logFC Pearson r (ON vs OFF) | 0.9932 | 去污染不扭曲表达估计 |
  | Tier A 保留率 (ON→OFF) | 96.77% | 真实信号大部分保留 |
  | Tier A Jaccard | 0.808 | ON/OFF 核心交集 4,587 基因 |
  | OFF-only 显著基因 | 1,038 | 去污染 OFF 多出的基因（潜在假阳性） |
  | ON-only 显著基因 | 187 | 去污染 ON 独有的基因（潜在假阴性/救回） |
  | 宿主 reads 保留率 | 65.23% | 非无菌昆虫样本的正常值 |

  - Tier 转移矩阵关键数字：
    - `4,587`：ON Tier A ∩ OFF Tier A（核心共识）
    - `657`：ON unclassified → OFF Tier A（最大转移块，污染假阳性主要来源）
    - `53`：ON Tier A → OFF unclassified（可能的假阴性）
  - 保持语言纪律：❌ ground truth，✅ reference-condition
  - 保持局限性声明

  **Must NOT do**:
  - 不要删除现有的 subsampling / degradation / ablation 段落
  - 不要声称 decontam 是 "validation"

  **QA Scenarios**:
  ```
  Scenario: 摘要包含所有 5 段落 + decontam 新数据
    Tool: bash
    Steps: grep 检查 benchmark_summary.md
    Expected: 含 "Tier A retention" 含 "96.77%" 含 "0.9932"
    Evidence: .sisyphus/evidence/benchmark-summary-check.txt
  ```

  **Commit**: NO (output artifact)

- [ ] 3. 执行 Task 8 一致性检查

  **What to do**:
  - `grep -i "ground.truth" results/benchmark/benchmark_summary.md` → 应无匹配
  - `grep -i "gold.standard" results/benchmark/benchmark_summary.md` → 应无匹配
  - 验证 5 个段落都存在（Annotation Degradation, Subsampling Stability, Decontamination, Ablation, Limitations）
  - 验证局限性声明包含：Drosophila 是模式生物、注释完整时共识优势缩小、需非模式昆虫验证
  - 验证 Bombyx case study 在摘要中被提及
  - 检查所有 Figure 引用路径指向实际存在的文件

  **Must NOT do**:
  - 不要修改 benchmark 脚本本身
  - 不要添加新的 benchmark 维度

  **QA Scenarios**:
  ```
  Scenario: 禁用词汇检查通过
    Tool: bash
    Steps: grep -ci "ground.truth\|gold.standard" results/benchmark/benchmark_summary.md
    Expected: 0（不区分大小写）
    Evidence: .sisyphus/evidence/benchmark-language-check.txt

  Scenario: 五段落全部存在
    Tool: bash
    Steps: grep -c "## " results/benchmark/benchmark_summary.md
    Expected: ≥ 5 个二级标题
    Evidence: .sisyphus/evidence/benchmark-sections-check.txt

  Scenario: Figure 引用路径有效
    Tool: bash
    Steps: 提取 benchmark_summary.md 中的 figures/ 路径，逐一检查文件存在
    Expected: 全部存在
    Evidence: .sisyphus/evidence/benchmark-figures-check.txt
  ```

  **Commit**: NO (validation pass)

---

## Execution Strategy

```
Task 1 (合并文件) ──┐
                    ├── Task 2 (更新摘要) ── Task 3 (一致性检查)
                    │
                    (所有文件操作，无依赖冲突)
```

三个 task 串行执行：先合并文件，再更新文档，最后检查。无并行性需求。

**Commit Strategy**: 无 commit（所有产出都是 output artifacts，按仓库规范不提交 `results/`）

---

## Success Criteria

### Verification Commands
```bash
# 1. 文件完整性
ls results/benchmark/*.tsv | wc -l     # 预期 ≥ 10
ls results/benchmark/figures/*.pdf | wc -l  # 预期 ≥ 11

# 2. 语言纪律
grep -ci "ground.truth\|gold.standard" results/benchmark/benchmark_summary.md
# 预期: 0

# 3. 关键数据存在
grep -c "96.77" results/benchmark/benchmark_summary.md   # 预期 ≥ 1
grep -c "0.9932" results/benchmark/benchmark_summary.md  # 预期 ≥ 1
grep -c "1038" results/benchmark/benchmark_summary.md    # 预期 ≥ 1
```

### Final Checklist
- [ ] `results/benchmark/` 包含全部 Task 1-6 输出
- [ ] `benchmark_summary.md` 包含 decontam COMPARE 数据
- [ ] 全文无 "ground truth" / "gold standard"
- [ ] 五段落完整
- [ ] 局限性声明包含 Bombyx case study 引用
- [ ] 所有 Figure 路径有效
