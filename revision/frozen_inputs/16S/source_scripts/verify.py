from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
import csv,json,collections,statistics,hashlib,zipfile
p=Path(__file__).resolve().parents[1]
def records(name):return {r.id:str(r.seq).upper() for r in SeqIO.parse(p/name,'fasta')}
def read(name):return list(csv.DictReader((p/name).open()))
x=read('accession_audit.csv');meta=read('species_panel.csv');allseq=records('all_eligible.fasta');rep=records('representatives.fasta');al=records('alignment.fasta');s=json.loads((p/'summary.json').read_text());assert set(allseq)=={r['accession'] for r in x if r['status']=='eligible'};assert set(rep)=={r['accession'] for r in meta};assert len({r['organism'] for r in meta})==len(rep);assert {a.removeprefix('_R_') for a in al}==set(rep);assert len({len(v) for v in al.values()})==1
orient=[]
for a,v in al.items():
 k=a.removeprefix('_R_');rev=a.startswith('_R_');assert v.replace('-','')==(str(Seq(rep[k]).reverse_complement()) if rev else rep[k]);orient.append((k,'reverse_complement' if rev else 'as_extracted'))
with (p/'alignment_orientation.csv').open('w') as f:
 w=csv.writer(f);w.writerow(['accession','alignment_orientation']);w.writerows(orient)
families=collections.Counter(r['family'] or 'unassigned' for r in meta)
with (p/'family_summary.csv').open('w') as f:
 w=csv.writer(f);w.writerow(['source_family','representative_species']);w.writerows(families.most_common())
with (p/'mitochondrial_records_to_review.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(x[0]));w.writeheader();w.writerows(r for r in x if r['mitochondrial_source_explicit']=='1' and r['status']!='eligible')
lens=list(map(len,rep.values()));s.update(representative_length_min=min(lens),representative_length_median=statistics.median(lens),representative_length_max=max(lens),alignment_columns=len(next(iter(al.values()))),reversed=sum(v=='reverse_complement' for _,v in orient));(p/'summary.json').write_text(json.dumps(s,indent=2))
(p/'README.md').write_text(f"""# 海葵目（Actiniaria）线粒体 16S 参考序列

NCBI Nucleotide 检索快照：2026-11-04。检索命中并下载 {s['downloaded']} 条 accession.version，保留 {s['eligible']} 条线粒体16S候选记录。具名物种面板含 {len(rep)} 个来源双名、{s['genera']} 属和 {s['families']} 个可识别来源科名。未定种记录可进入全部合格序列库，但不计入物种代表。分类使用 GenBank 来源标签，未对 WoRMS 同物异名进行合并。

- [全部合格序列](all_eligible.fasta) · [每种一条代表](representatives.fasta) · [物种与 accession 映射](species_panel.csv) · [各科统计](family_summary.csv)
- [代表序列比对](alignment.fasta) · [比对方向](alignment_orientation.csv) · [逐 accession 审核](accession_audit.csv) · [有线粒体来源信息的待复核记录](mitochondrial_records_to_review.csv) · [检索式、UID与原始分批 GenBank](downloads/)

提取优先采用16S/rrnL的rRNA或gene注释，按feature坐标、链方向与拼接结构提取。对仅标记为大亚基的注释要求线粒体来源证据；对明确16S的单基因短记录允许标题回退。筛选要求至少120 bp、模糊碱基≤5%。{s['method_counts'].get('feature_rRNA',0)}条来自rRNA注释，{s['method_counts'].get('single_locus_title',0)}条依据单基因标题。无可识别16S边界的记录保留原始GenBank与排除原因，未将整条线粒体基因组当作16S序列。

每个具名来源物种只保留一条。优先采用完整基因组中明确注释边界、且覆盖 Metridium senile 线粒体参考 NC_000933.1 至少90%明确16S位点的记录。其次按比对至该参考的明确碱基位点数、序列自身明确碱基数、较少模糊碱基、注释来源、RefSeq标记及accession依次排序。所有候选及其同种排名见[筛选明细](selection_ranking.csv)。复现顺序为 curate.py、select_complete.py、verify.py。代表片段范围 {min(lens)}–{max(lens)} bp，中位数 {statistics.median(lens):g} bp。MAFFT 使用 --adjustdirection --retree 2 --maxiterate 0；比对长 {s['alignment_columns']} 列。_R_ 前缀标记相对于提取序列反向互补的记录。已逐条验证去除比对空位后可恢复对应方向的原序列。部分代表仅覆盖16S的一段，后续引物设计应按每个结合窗口重新确定可评价物种。
""")
files=sorted(f for f in p.rglob('*') if f.is_file() and f.name!='SHA256SUMS.json' and '__pycache__' not in f.parts);(p/'SHA256SUMS.json').write_text(json.dumps({str(f.relative_to(p)):hashlib.sha256(f.read_bytes()).hexdigest() for f in files},indent=2))
with zipfile.ZipFile(p.parent/'Actiniaria_16S.zip','w',zipfile.ZIP_DEFLATED) as z:
 for f in files+[p/'SHA256SUMS.json']:z.write(f,str(f.relative_to(p.parent)))
with zipfile.ZipFile(p.parent/'Actiniaria_16S.zip') as z:assert z.testzip() is None
print(s)
