# 已知问题（待后续修复）

**生成日期**: 2026-05-05
**分支**: discussion/consensus-downstream

---

## 高优先级

### 无测试基础设施
- 零 test 文件，零 CI/CD
- 验证全靠手动 `snakemake --dry-run`
- 建议: 添加 `pytest`（Python 脚本）+ `testthat`（R 脚本），至少覆盖数据转换函数

### 无容器化
- 无 Dockerfile / docker-compose.yml
- 可重复性依赖裸机 Conda 环境
- 建议: 编写 Dockerfile，或在 README 中提供 `conda-lock` 生成指令

### 无版本管理
- 无 git tag，无 changelog，无 `version.py`
- 建议: 对稳定版本打 tag（如 `v1.0.0`），维护 `CHANGELOG.md`

### 无 LICENSE
- 无许可证文件
- 法律上他人无法合法复现或修改

---

## 中优先级

### 仅中文文档
- README、注释、commit message 全中文
- 限制国际合作和社区贡献
- 建议: README 提供英文版本，关键技术注释双语标注

### 无示例数据
- `data/` 目录为空，新用户无法快速跑通验证
- 建议: 捆绑小型测试数据集（如 2 样本 × 2 条件，子采样 reads，迷你参考基因组）

### .snakemake/ 运行时目录泄漏
- Snakemake 运行时产物进入工作目录
- 虽然 .gitignore 排除了，但不干净
- 建议: 配置 `snakemake --directory` 或在文档中说明

### qc.yaml 通道数不一致
- `qc.yaml` 有 3 个通道（含 `pkgs/r`），其他 8 个文件只有 2 个
- 功能上无影响（R 包已通过 conda-forge/bioconda 覆盖），但形式上不统一

---

## 低优先级

### 无 linter/formatter 配置
- 无 black/ruff（Python），无 lintr（R），无 pre-commit
- 建议: 添加 `.pre-commit-config.yaml` + `pyproject.toml` 配置 ruff

### Conda 环境未锁定
- 无 `conda-lock.yml`，依赖随时间漂移
- 建议: 每次环境变更后生成 `envs/*.lock.yaml`

### DESeq2 锁定，无 edgeR/limma 切换入口
- 低生物学重复数场景下 edgeR 可能更合适（DESeq2 对低重复数较保守）
- 建议: 评估是否需要为低重复数场景保留 edgeR 入口

---

## 已修复（历史反模式记录）

以下问题曾在 `run_consensus_dea.R` 中存在，已于 2026-05-05 修复。保留此记录以防止回退。

### `choose_direction()` 用 P 值断案（胜者诅咒）

**旧逻辑**（L184-189 原代码）：
```r
if (rra_up_p < rra_down_p) return("up")
if (rra_down_p < rra_up_p) return("down")
```
**问题**：上调/下调支持数平局时，选择 P 值更低的方向。被选中的低 P 值随后被用作 `best_rra_fdr` — 循环依赖。等价于两次检验取 min。

**修复**：删除 P 值断案。平局 → `mixed`。同时引入 ≥2× 优势规则，边际冲突（如 2 up vs 1 down）也判 `mixed`。

### Mixed 基因用 `pmin()` 取最小 P 值

**旧逻辑**（L617-620 原代码）：
```r
best_rra_fdr[consensus_direction == "mixed"] <- pmin(rra_up_fdr, rra_down_fdr)
```
**问题**：对两个相反假设的 P 值取最小值，人为制造假阳性。虽然这些基因不进入 Tier（tier == "Conflict"），但 `consensus_results.tsv` 中的 p 值仍有误导。

**修复**：mixed/none 方向基因的共识 P 值设为 `NA`。删除四行 `pmin()`。

### Mixed 基因作为单独 "Conflict" Tier

**旧逻辑**：`consensus_direction == "mixed"` → `tier = "Conflict"`
**问题**："Conflict" 暗示这些基因可能仍有部分证据。实际上它们连方向都没法定，不应该在 tier 中出现。

**修复**：→ `tier = "unclassified"`，与 support_n < 2 的基因同样处理。

### `exact = length(rank_lists) < 10` 动态切换精确模式

**旧逻辑**：
```r
exact = length(rank_lists) < 10  # 4 定量器 → exact = TRUE
```
**问题**：`?rhoScores` 文档明确警告：*"exact: whether exact p-value should be calculated. Warning: computationally unstable and does not give considerable gain."* 精确模式内部使用 Stuart 方法组合 — 而 Stuart 方法原作者明确说"不能用于显著性判定"。

**修复**：→ `exact = FALSE`，始终使用 Bonferroni 校正。

### "RRA 主证据 / CCT 次证据" 的措辞

**问题**：CCT 在数学上对相关 P 值有完备保证，RRA 反而有独立性假设缺陷。称 CCT 为"次证据"在审稿人面前站不住脚。

**修复**：→ "RRA 与 CCT 作为并列的双重共识引擎（co-primary）。进入 Tier A 必须两者同时满足。" 代码注释和 `docs/methods.md` 同步更新。

### consensus_summary.tsv 缺少 RRA 方法敏感性指标

**问题**：无法回答审稿人"结果对 rank aggregation 方法选择敏感吗？"的疑问。

**修复**：新增 `rra_mean_tier_a/b/c_concordance` 指标（基于 Borda mean 聚合的替代 Tier 判定），写入 `consensus_summary.tsv` 和 `sensitivity_analysis.tsv`。
