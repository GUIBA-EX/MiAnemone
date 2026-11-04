from pathlib import Path
import csv,re
ROOT=Path(__file__).resolve().parents[1]
def rows(p):return list(csv.DictReader((ROOT/p).open()))
body=(ROOT/'publications/MiAnemone_preprint.md').read_text();banner='![MiAnemone overview](figures/MiAnemone_overview.png)\n\n';en=body.replace('](figures/','](publications/figures/').replace('](MiAnemone_supplement.md)','](publications/MiAnemone_supplement.md)');en=en.replace('## Abstract','[中文](README.md) · [Revision PDF](publications/MiAnemone_preprint.pdf) · [Response to review](publications/response_to_review.md) · [Reproduction guide](revision/README.md)\n\n## Abstract',1);(ROOT/'README.en.md').write_text(banner+en)
zh='''![MiAnemone 概览](figures/MiAnemone_overview.png)

# MiAnemone: in silico design and evaluation of mitochondrial 12S and 16S rRNA primers for sea anemone (Actiniaria) eDNA metabarcoding

**Fei Xia** · 2026-11-04 计算预印本 · 第二轮修订

[English](README.en.md) · [修订论文 PDF](publications/MiAnemone_preprint.pdf) · [逐条回复](publications/response_to_review.md) · [复现说明](revision/README.md)

## 摘要

本研究对海葵目线粒体12S和16S rRNA的22对普通与简并引物进行计算评估。原设计面板分别包含250和265个具名物种。本轮补充公共物种面板、非嵌套普通/D的因子分解、错配位置、统一目标/背景模型、全部具名记录敏感性、固定窗口属级留出及条形码区分度。12S XS的G10A使完全匹配从2增至192/245，覆盖仍为205种，但该替换使用了事后放宽的温差门槛。原候选中三组D还改变了固定碱基。现已用嵌套候选替换12S XSD，覆盖227/245、完全匹配212/245，严格背景仍为9条。现有结果支持不同用途下的计算候选优先级，本研究建立了后续试剂评价的计算优先级。

## 引言

海葵目的分类研究曾报告约1,200个描述种，并使用12S和16S等标记重建亲缘关系（Rodríguez等，2014）。本研究借鉴MiFish的保守侧翼思路和Anth-28S的简并策略，分别考察结合位点兼容性与内部条形码的信息量。珊瑚虫线粒体演化较慢（Shearer等，2002；Hellberg，2006），较高覆盖并不自动带来较高物种区分度。引物错配可影响扩增偏差（Piñol等，2015），但其作用还取决于位置与反应条件，因此完全匹配增加不能直接换算为检出率增加。

## 材料与方法

### 参考库与比对

原始检索分别获得1,360条12S和2,527条16S记录，可用记录池为786和1,245条。采用GenBank物种名称与谱系，排除sp./cf.等非明确双名记录后，每种选择一个参考，形成250种12S和265种16S设计面板。12S实际排名按明确碱基数、rRNA注释、RefSeq和登录号；12S排序未使用基因组完整性。16S另使用基因组完整性与参考列覆盖排名。尚未完成WoRMS同物异名全面整理或标本误鉴定核验。

MAFFT 7.526参数为`--adjustdirection --retree 2 --maxiterate 0 --thread 4`，不是结构感知比对。坐标为*Metridium senile* NC_000933.1相应基因的一起算闭区间；12S位于线粒体12,201–13,282，16S位于72–2,260。普通/D共享四个端点，但不保证D必然包含普通序列。

### 候选生成、热力学与替换规则

初始搜索每3个参考碱基取一个起点，尝试20/22 nt，要求至少10种或10%目标物种提供完整位点。普通候选来自观察模体与共识；简并位点限于3′半段，每端至多8个变体。按覆盖、完全匹配、较低简并度排序，并在坐标分箱内保留有限候选，再于所选窗口重新配对。该过程有候选截断，是启发式搜索而非全局最优证明。

末端过滤排除可形成五连同碱基、末三位同碱基及不符合GC条件的候选。末位可简并，但兼容须由实际变体提供。primer3-py 2.4.0使用50 mM单价盐、1.5 mM二价盐、0.6 mM dNTP，分别评估50 nM/变体和200 nM/方向均分两种情景。初始最低Tm为45°C、结构Tm低于55°C、F/R温差≤5°C，没有额外全序列GC百分比门槛。

G10A所用6°C例外是检查该候选后提出的事后规则。补充表列出了23个后续候选在ΔTm 4–7°C及新碱基频率80–95%网格下的全部判定。正式22对序列未因本轮解释性控制而自动替换。

### 配对模型与评价指标

可评价要求参考跨越两端且结合区没有模糊碱基，允许比对缺口；完整等长模体是更严格子集。原始双链扫描要求产物100–500 bp、内部≥50 bp；覆盖为每端≤2错配且末三位零错配，完全匹配为两端零错配。D允许任意一个F变体与任意一个R变体共同匹配。

原背景为322个12S刺胞动物物种和1,775条16S刺胞动物记录，原扫描允许内部0–1,000 bp。严格WC/GT使用两错配及末三位规则；宽松WC为三错配无末端限制，宽松GT为两错配无末端限制。GT将指定的潜在摇摆配对计为兼容，包括末端。本轮额外统一目标和背景几何范围，比较WC/GT、每端0–4错配、末0/3/5位保护。位置连续罚分是自定义诊断，不是校准的PCR效率。

公共覆盖面板取各窗口均可评价的物种；公共区分度进一步要求各窗口内部序列全部明确。计算物种唯一性、覆盖∩唯一性、最近邻编辑距离、属/科归属、多记录观察barcode gap和留一记录最近邻分类。Wilson区间以当前面板为条件，不代表海葵目随机抽样。

全部具名记录池为694条12S和1,116条16S，报告任一/多数/全部记录匹配及种内位点变异。五折按属分块，在训练属内重建普通/D并经过热力学过滤，再扫描留出属。端点与比对仍来自全数据，背景未用于调节这些诊断候选，因此属于固定窗口条件测试，不是完全独立的窗口发现验证。

## 结果

### 公共面板与非嵌套对照

公共覆盖面板为12S的240种和16S的67种。16S SD、LD、MD和XLD分别覆盖62、59、58和58/67种；HD为65/67，但其宽松背景较高，不应仅按覆盖推荐。各窗口原始分母不同，不能直接将其百分比当成同一组物种排名。

![公共面板覆盖与完全匹配](publications/figures/review_common_panel.svg)

**图1.** 各标记共同可评价物种面板中的覆盖与完全匹配。柱上为物种数，误差线为条件Wilson 95%区间；H/HD作为历史位点对照置于最后。

![参考位点变异、有效物种数与引物位置](publications/figures/review_entropy.svg)

**图2.** 每个标记依次展示位点Shannon熵、贡献明确A/C/G/T的物种数、分行排列的F/R结合区间，三层共享NC_000933.1对应基因的一起算坐标。缺口和模糊碱基不计入分母；少于两个有效物种的位置不绘制熵值。低样本量下的低熵不表示已证实保守。箭头覆盖实际结合区间并指向3′端，连接线表示中间区域。普通与D端点相同，共用一行；12S与16S使用各自的坐标尺度。

替换前的三组非嵌套对是12S XS/XSD、12S S/SD和16S XS/XSD。仅将简并位点加入普通12S XS，覆盖为227/245、完全匹配为212/245；仅改变固定碱基则为205覆盖、0完全匹配。两类变化同时存在才得到原XSD的230覆盖、0完全匹配。因此零完全匹配主要是固定碱基效应。16S XS嵌套控制为223覆盖、167完全匹配，但末端过滤不通过。12S XS嵌套对照两项过滤均通过，最大温差4.84°C，现正式替换为XSD并列为短片段优先候选，普通XS保持不变。

### 共识替换与位置效应

12S XS的230个完整模体中，G10位点有228个A、2个G；另15个可评价物种存在缺口或结合段长度变化。G10A后的F/R温差为50 nM时5.18°C、200 nM时5.33°C。12S S的F9C→A保持173种覆盖，完全匹配0→4，最大温差5.95°C，未满足事后规则的10个百分点增益。L的F1T→C使完全匹配1→202，但两种Lebrunia不再符合两错配规则；原XSD另涉及Phyllodiscus，均属归档分类中的Aliciidae。这是阈值预测，不是必然不能PCR扩增。

![错配距3′端的位置分布](publications/figures/review_mismatch_profiles.svg)

**图3.** 明确同源碱基中的错配比例，横轴从3′末端的0位开始；竖线分隔末端三位。分母排除缺口及模糊碱基，白格表示超出引物长度或没有明确观测。

### 统一背景模型、种内记录与留出

统一长度后，16S H/HD的宽松WC背景从486/976变为452/939，仍然较高；LD由46降至45。原12S XS、XSD、M、MD、L和LD的严格WC命中是完全相同的9条黑珊瑚记录。两套扫描实现在统一严格WC/GT模型下完成57,584次记录—模型比较，零差异。原500/501计数没有截断。

![统一模型下目标与背景变化](publications/figures/review_unified_models.svg)

**图4.** 相同产物几何约束下的目标覆盖与背景命中数，末三位保持匹配，每条曲线连接每端允许0至4个错配的结果。横轴采用非负symlog尺度，0–1为线性区；实线为普通、虚线为D。

![全部引物的分类群背景命中](publications/figures/review_background_all.svg)

**图5.** 22对引物在原始背景面板中的命中数，按分类群分组。W/G为严格WC/GT，W+/G+为宽松模型；仅标注非零值，两幅面板共享同一计数色阶。各标记背景库大小不同，原始计数不能直接解释为标记间特异性差异。

全部记录分析中，12S XS在247个可评价物种里，任一/多数/全部记录覆盖为208/206/203；16S LD为233/232/225，分母259。重新比对也影响可评价分母。属级重选中，12S M普通/D分别覆盖219/226与219/226、完全匹配197与205；16S L普通/D分别覆盖137/246与213/246、完全匹配106与186。候选在每折重新生成，不能解释为对固定正式引物的外部验证。

### 区分度、结合位置与使用限制

原12S唯一内部序列为XS 35/244、S 44/244、M 86/239、L 68/243、XL 100/239。L低于M，说明区分度不随长度单调增加。共同完整面板进一步将扩增兼容性与唯一性相交：

![共同面板覆盖且唯一](publications/figures/review_joint_resolution.svg)

**图6.** 共同完整插入序列面板中，同时满足覆盖与参考序列唯一性的物种比例，柱上为物种数。该指标描述已采样参考库，不等同于实验验证的物种鉴定成功率。

多记录分析中，12S XS只有11/111个多记录物种具有正的观察barcode gap，XL为28/110；严格留一记录最近邻正确归属为52/658和179/646。这些结果仍受种内采样不足影响。

12S XS-F与M/MD-R、SD-F与M/MD-R等共享区域的互补作用使跨配对异二聚体Tm超过55°C，不宜直接混为多重体系。实际产物的最小/中位/最大长度和名义双端20 bp重叠统计见论文Table2及`revision/length_summary.csv`。实际接头待定，尚未评估带接头结构。

## 使用优先级与讨论

| 用途 | 计算优先候选 | 主要取舍 |
| --- | --- | --- |
| 短片段检测 | 12S XSD（采用Table3，0/1） | 末端及热力学过滤通过，覆盖227/245、完全匹配212/245；未评价鱼类、甲壳类共扩增 |
| 覆盖与区分度兼顾 | 12S MD | 完全匹配较高；存在黑珊瑚背景 |
| 较广16S目标兼容性 | 16S LD | 较L显著增加覆盖；宽松背景也增加 |
| 黑珊瑚共存环境 | 对本地参考比较12S S/SD及16S候选 | 背景与目标完整性必须一起判断 |
| 手工位点历史比较 | H/HD | 宽松背景过高，不优先推荐 |
| 多重PCR | 不作现成混合推荐 | 共享位点二聚体及接头需独立评估 |

同坐标隔离的是位置差异，只有嵌套序列才能进一步解释简并效应。覆盖、完全匹配和物种分辨率需要分别报告。当前模型不是环境样本检出率；没有新增湿实验。背景结果限定于原刺胞动物参考库。独立窗口发现、冻结后时间外推及分类名称全面核验尚未完成。

## 引物与逐碱基匹配图

正式序列见[FASTA](primers.fasta)；每条oligo的绝对Tm、GC、长度与简并度见[属性表](revision/oligo_properties.csv)。原有逐碱基匹配图仍保留，下列每一对均可单独查看。

'''
# Existing diagrams are linked by actual filenames only.
paths=sorted((ROOT/'figures').glob('*'))
for gene in ['12S','16S']:
 hits=[p for p in paths if gene in p.name and ('base' in p.name.lower() or 'match' in p.name.lower())]
 if hits:zh+='\n'+gene+'：'+' · '.join(f'[{p.stem}](figures/{p.name})' for p in hits)+'\n'
zh+='\n'+' · '.join(f'[{q.stem}](figures/template_matching/{q.name})' for q in sorted((ROOT/'figures/template_matching').glob('*.svg')) if '_preview' not in q.name)+'\n'
zh+='\n| 基因 | 配对 | F（5′–3′） | R（5′–3′） |\n| --- | --- | --- | --- |\n'
for r in rows('all_results.csv'):zh+='| '+' | '.join(r[k] for k in ['gene','name','F','R'])+' |\n'
refs=body.split('## References\n',1)[1].split('## Table 1.',1)[0]
zh+='\n## 数据与复现\n\n[GitHub项目](https://github.com/GUIBA-EX/MiAnemone) · [分析索引与复现命令](revision/README.md) · [逐条审稿回复](publications/response_to_review.md)。本目录是独立工作稿，没有将其标记为已公开的新版本。作者单位、ORCID、资助及利益冲突声明均待定，未代填。著作权声明沿用[COPYRIGHT.md](COPYRIGHT.md)。\n\n## 参考文献\n\n'+refs
(ROOT/'README.md').write_text(zh)

for file in ['README.md','README.en.md']:
 p=ROOT/file;s=p.read_text();lines=s.splitlines();idx=next(i for i,l in enumerate(lines) if l.startswith('**Fei Xia**'))
 lines[idx]='**Fei Xia**';lines[idx+1:idx+1]=['','<!-- publication:start -->','<!-- publication:end -->',''];s='\n'.join(lines)+'\n'
 s+='\n## '+('引用 / Citation' if file=='README.md' else 'Citation')+'\n\nXia, F. (2026). Computational preprint.\n'
 s+='\n'+('[著作权声明](COPYRIGHT.md)' if file=='README.md' else '[Copyright and permissions](COPYRIGHT.md)')+'\n';p.write_text(s)
