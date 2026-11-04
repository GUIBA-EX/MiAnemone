# MiAnemone 审稿修订分析

本轮修订日期：2026-11-04。12S XSD已采用通过末端与热力学过滤的嵌套候选，普通XS不变。当前统计与图表均已按新序列更新；prior_XSD保存替换前的因子对照和序列。

## 数据与分析索引

| 文件 | 用途与边界 |
| --- | --- |
| panel_comparison.csv | 原窗口与公共物种面板的覆盖、完全匹配和覆盖率 Wilson 区间；各目不是随机抽样 |
| joint_resolution.csv | 公共完整插入序列面板上的唯一性与覆盖∩唯一性 |
| nesting.csv / factorial_controls.csv | 当前普通/D的嵌套关系及两组非嵌套分解；prior_XSD/factorial_filters.csv保留Table3的三组替换前对照 |
| mismatch_positions.csv / site_sensitivity.csv | 距 3′ 端的错配位置及阈值敏感性；位置罚分是探索性指标 |
| family_coverage.csv / template_status.csv | 按科覆盖与每条模板的可评价、缺口及长度变化 |
| threshold_sensitivity.csv | 23 个候选在 ΔTm 4–7°C、新碱基频率 80–95% 下的判定；其他门槛固定 |
| oligo_properties.csv | 每条 oligo 的长度、GC、变体数、两种浓度情景下绝对 Tm |
| cross_assay_dimers.csv | 同一基因不同配对之间的变体异二聚体；尚无实际接头 |
| unified_models.csv / unified_record_hits.csv | 目标与原背景使用相同长度约束；WC/GT、每端 0–4 错配、末 0/3/5 位保护 |
| independent_validation.json | 两套实现对统一严格模型的独立检查 |
| product_lengths.csv / length_summary.csv | 实际预测产物及 20 bp 名义双端重叠可行性；尚未扣除质量剪切 |
| nearest_neighbors.csv | 编辑距离最近邻、属/科一致性；非 BLAST/IDTAXA 分类器 |
| all_record_summary.csv / all_record_species.csv | 全部具名记录的任一/多数/全部记录覆盖与位点变异；重新比对也影响分母 |
| observed_barcode_gap.csv | 多记录物种的观察到的种内最大/种间最小距离，以及留一记录最近邻归属 |
| redesign_candidates.csv / redesign_holdout.csv | 五折属级条件留出；训练内重建普通/D，固定已选端点，未重新发现窗口，未用背景调参 |
| genus_holdout.csv | 早期仅按观察到的完整配对模体重选的诊断，不作为正文热力学重选结果 |
| frozen_inputs/ | 下载归档、原始筛选记录和源脚本快照 |
| environment.json | 实际运行环境版本 |

## 复现

在仓库根目录运行。依赖版本见 `environment.json`；MAFFT、C++ 编译器、pandoc 和 Chrome/Chromium 需可执行。Python 依赖为 Biopython、numpy、primer3-py、matplotlib。

```sh
python3 scripts/review_analysis.py
python3 scripts/review_models.py
python3 scripts/review_extensions.py
python3 scripts/review_cross_validation.py
python3 scripts/review_all_records.py
c++ -O3 -std=c++17 scripts/review_exhaustive.cpp -o /tmp/mianemone_review_scan
/tmp/mianemone_review_scan revision/independent_12S
/tmp/mianemone_review_scan revision/independent_16S
/tmp/mianemone_review_scan revision/all_records_12S
/tmp/mianemone_review_scan revision/all_records_16S
python3 scripts/review_record_summary.py
python3 scripts/review_barcode_gap.py
python3 scripts/check_factorial_filters.py
python3 scripts/check_revision_consistency.py
python3 publications/plot_revision.py
python3 publications/assemble_revision.py
python3 publications/sync_revision_readmes.py
python3 scripts/stamp_preprint.py 1.4
python3 publications/build_preprint.py
```

原始筛选脚本快照保留其历史路径假定，用于追溯生成过程；上面的修订分析使用当前仓库内数据。`review_exhaustive.cpp` 仅将 WC 严格目标产物写入条形码表，GT 命中单独保存在 hit 表。`independent_validation.json` 的模型范围为 WC/GT、每端≤2错配、末三位保护及相同产物几何范围。

## 本轮未完成或暂缓项目

- 冻结日期之后的时间外推：没有声称获得未来序列测试集。
- 独立窗口发现：已有固定窗口的属级重选；尚未完成训练内重新选择全部窗口的嵌套验证。
- WoRMS 同物异名全面核对及标本误鉴定排查：仍沿用归档 GenBank 名称。
- 实际接头结构：接头序列待定。
- 作者单位、ORCID、资助、利益冲突声明：作者回复“待定”，未填造。
- DOI：尚未申请；本次结果由GitHub预印本版本标签固定。

这些项目不被标记为已解决，逐条回复见 `../publications/response_to_review.md`。
