# Prochl-SBE 仓库用途与 BadiRate / AnGST / GeneRax 上手指南

## 1. 这个仓库在做什么？

这个仓库主要是 **分子定年分析（molecular dating）及其上游数据准备**，核心包括：

1. `HomoTest/`：筛选组成成分更均一（compositionally homogeneous）的直系同源基因家族。
2. `phylogeny_construction/`：对同源基因家族做比对、裁剪、拼接与分区建模（MAFFT、trimAl、PartitionFinder）。
3. `dRdC_cal/`：计算 dr/dc 相关指标。
4. `mcmctree_without_outgroup/` 与 `Mcmctree_with_outgroup/`：给出 MCMCTREE 三步流程示例（BASEML 估速率 -> MCMCTREE 估 Hessian/gradient -> MCMCTREE 跑后验时间）。

> 注意：仓库本身没有直接提供 BadiRate、AnGST、GeneRax 的脚本，但提供了这些方法所需的大部分基础输入（基因家族比对、拼接树构建流程、时间树流程示例）。

---

## 2. 先把仓库原有流程跑通（建议）

以下命令用于获得可用于后续三种方法的标准化输入。

### 2.1 基因家族比对与清洗（在 `phylogeny_construction/`）

```bash
cd phylogeny_construction
perl stp1.mafft.pl
perl stp2.format.align.pl
perl stp3.trimAl.pl
perl stp3.format.trimAl.pl
```

### 2.2 生成拼接矩阵与分区文件

```bash
perl stp5.generate_data.for_cat.pl
perl stp6.partition_finder_pre.a.pl .
bash stp7.partition_finder.a.sh
```

执行后常用关键文件：
- `ptt_fdr.a/concat.phy`（拼接比对）
- `ptt_fdr.a/partition_finder.cfg`（分区定义）

### 2.3 复用仓库中的时间树示例（可选）

如果你要让 GeneRax/AnGST 在“已定年的物种树”上做分析，可参考：

```bash
cd ../mcmctree_without_outgroup/stp1.substitution_rate
bash baseml.sh

cd ../stp2.matrix
bash mcmctree.stp2.sh

cd ../stp3.mcmctree
bash mcmctree.stp3.sh
```

---

## 3. 三种方法怎么接入

下面给出“如何把仓库数据接入三种工具”的统一思路。

## 3.1 GeneRax（基因树-物种树联合 + DTL）

GeneRax 典型输入是：
- 物种树（Newick）
- 每个基因家族的 MSA（FASTA/PHYLIP）
- 基因-物种映射（mapping）

### 推荐做法

1. 用 `phylogeny_construction/` 产出的 `*.mnf.trimAl` 作为每个 family 的 MSA。
2. 物种树可用你已有的参考树（例如 MCMCTREE 使用的树拓扑）。
3. 基因 ID 命名统一成 `species|gene`（仓库脚本中已有类似格式处理），写成 mapping 文件。

### 示例命令（按你本机 GeneRax 版本调整）

```bash
generax \
  --species-tree species_tree.nwk \
  --families families.txt \
  --rec-model UndatedDTL \
  --prefix generax_out \
  --strategy SPR
```

其中 `families.txt` 一般每个 family 指向：
- 对应 MSA 路径
- 对应初始基因树（可由 IQ-TREE/RAxML 预估）
- 映射文件

---

## 3.2 AnGST（基因树与物种树 reconciliation）

AnGST 通常需要：
- 物种树
- 每个基因家族树（可带 bootstrap 样本）
- 事件代价设置（dup/loss/transfer）

### 推荐做法

1. 基于 `*.mnf.trimAl` 给每个家族建基因树（例如 IQ-TREE）。
2. 准备一个全局物种树（可用 MCMCTREE 输入树或其它参考树）。
3. 在 AnGST 配置里设定事件代价（可先用默认，再做敏感性分析）。

### 示例命令（示意）

```bash
python angrst.py \
  --species species_tree.nwk \
  --genes gene_trees/ \
  --costs "D=2,L=1,T=3" \
  --output angrst_out
```

> 注：AnGST 各发行版参数名可能不同，请以你安装版本文档为准。

---

## 3.3 BadiRate（基因家族获得/丢失速率）

BadiRate 一般关心：
- 家族拷贝数矩阵（species × family）
- 物种树（可超度量树或分支长度树，视模型设定）

### 推荐做法

1. 从同源家族结果统计每个物种在每个 family 的拷贝数，形成 `count_matrix.tsv`。
2. 使用统一物种树（和上面 GeneRax/AnGST 一致，便于比较）。
3. 先跑一个全局速率模型，再跑分支/谱系异质模型比较 AIC/BIC。

### `count_matrix.tsv` 具体怎么得到？

可直接用仓库里的家族比对文件（例如 `phylogeny_construction/*.mnf.trimAl`）按序列头统计拷贝数。
本仓库新增了一个脚本：`tools/build_count_matrix.py`。

> 统计规则：若序列头是 `>species|gene`，则按 `|` 左边作为物种名；同一个 family 里同一物种出现几条序列，就记为几拷贝。

```bash
# 在仓库根目录执行
python tools/build_count_matrix.py \
  --indir phylogeny_construction \
  --pattern "*.mnf.trimAl" \
  --out count_matrix.tsv
# 或者（仓库根目录下）
python build_count_matrix.py --indir phylogeny_construction --pattern "*.mnf.trimAl" --out count_matrix.tsv
```

如果你想固定物种顺序（并补齐缺失物种为 0），可提供物种列表（每行一个物种名）：

```bash
python tools/build_count_matrix.py \
  --indir phylogeny_construction \
  --pattern "*.mnf.trimAl" \
  --species-list phylogeny_construction/genome.list \
  --out count_matrix.tsv
```

输出文件结构示例：

```text
Species	COG0001	COG0002	...
Prochl_A	1	0	...
Prochl_B	2	1	...
```

拿到这个 `count_matrix.tsv` 后，就可以直接作为 BadiRate 的输入矩阵。

### 示例命令（示意）

```bash
badiRate \
  --tree species_tree.nwk \
  --counts count_matrix.tsv \
  --model birth-death \
  --out badirate_out
```

> 注：BadiRate 不同版本语法差异较大，上述为“参数骨架”。

---

## 4. 三方法结果如何统一比较（建议）

为了让三种方法可比，建议统一：

1. **物种树版本**：同一棵拓扑（最好同一套分支长度）。
2. **基因家族集合**：三种方法使用同一批 family。
3. **命名规则**：species/gene ID 一次性规范。
4. **事件定义对照表**：
   - GeneRax: D/T/L
   - AnGST: D/T/L
   - BadiRate: gain/loss（不直接给转移时，可与前两者联合解释）

最终可以做一张表：每个分支的 gain/loss/transfer 支持度和模型证据（AIC/BIC/likelihood）。

---

## 5. 一套可直接执行的最小工作流（模板）

```bash
# 0) 在 phylogeny_construction 中生成清洗后的家族比对
cd phylogeny_construction
perl stp1.mafft.pl
perl stp2.format.align.pl
perl stp3.trimAl.pl
perl stp3.format.trimAl.pl

# 1) 给每个 *.mnf.trimAl 建基因树（示例用 IQ-TREE）
mkdir -p ../gene_trees
for f in *.mnf.trimAl; do
  base="${f%.mnf.trimAl}"
  iqtree2 -s "$f" -m MFP -bb 1000 -nt AUTO -pre "../gene_trees/${base}"
done

# 2) 准备 species_tree.nwk + mapping 文件 + family 列表
#    (mapping 可按 species|gene 命名规则自动生成)

# 3) 跑 GeneRax / AnGST / BadiRate
#    命令按你的软件版本替换
```

如果你愿意，我可以下一步直接帮你：
- 基于你当前文件命名，自动生成 GeneRax `families.txt` 与 mapping；
- 从现有基因家族 fasta 自动统计 BadiRate 的拷贝数矩阵；
- 生成一份可复现的 `run_all.sh`（三工具顺序跑 + 汇总结果）。
