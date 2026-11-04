from choose import p,d,score,panel,cols,al,write
from independent import flags
from Bio import SeqIO
import collections,csv,json,hashlib,zipfile
v=p/'verify';rs=list(csv.DictReader((p/'comparison.csv').open()));ids={q['id']:q for q in rs};orig=rs[0]
cpp=collections.defaultdict(dict)
with (v/'independent_hits.tsv').open() as h:
 for x in csv.DictReader(h,delimiter='\t'):
  if x['name'] in ids:cpp[x['name']][x['role'],x['model'],x['accession']]=x
bg={r['accession']:r for r in csv.DictReader((d/'background_panel.csv').open())};tax={r['accession']:r for r in panel}
targets={r.id:str(r.seq).upper() for r in SeqIO.parse(d/'target.fasta','fasta')};back={r.id:str(r.seq).upper() for r in SeqIO.parse(d/'background.fasta','fasta')};checks=0;out=[];detail=[];fa=[];diag=[]
for q in rs:
 k=q['id'];name='H' if k=='H' else q['group']+('D' if max(int(q['deg_F']),int(q['deg_R']))>1 else '');q['name']=name
 coords=[cols[int(q[x])-1] for x in ('F_start','F_end','R_start','R_end')];a,b,c,e=coords
 counts=collections.Counter()
 for role,lib in [('target',targets),('background',back)]:
  for acc,s in lib.items():
   ok=True
   if role=='target':
    t=al[acc];ok=bool(t[:a+1].strip('-') and t[e:].strip('-') and set(t[a:b+1]+t[c:e+1])<=set('ACGT-'));counts['evaluable']+=ok
   for model in (['WC'] if role=='target' else ['WC','GT']):
    ex,cv,broad=flags(s,q['F'],q['R'],role=='target',model=='GT');h=cpp[k].get((role,model,acc),{})
    assert [ex,cv,broad]==[h.get(y)=='1' for y in ['exact','le2last3','le3' if role=='background' and model=='WC' else 'le2']],(k,role,acc,model)
    checks+=1
    if role=='target':
     counts['exact']+=ok and ex;counts['covered']+=ok and cv
     detail.append(dict(name=name,accession=acc,species=tax[acc]['organism'],family=tax[acc]['family'],evaluable=int(ok),covered=int(ok and cv),exact=int(ok and ex)))
    else:counts[model]+=cv;counts['WC3' if model=='WC' else 'GT2']+=broad
 for y in ['covered','exact','evaluable','WC','GT','WC3','GT2']:assert counts[y]==int(q[y]),(k,y)
 if k!='H':
  dest=diag if q['group']=='fixed' else fa
  dest.extend(f">{name}_{side}\n{q[side]}\n" for side in 'FR')
 out.append(q)
write(p/'candidate_results.csv',out);write(p/'species_results.csv',detail)
(p/'primers.fasta').write_text(''.join(fa));(p/'original_site_repair.fasta').write_text(''.join(diag))
(p/'validation.json').write_text(json.dumps(dict(record_model_checks=checks,disagreements=0,background_records=len(back),target_species=len(targets)),indent=2))
lines=['# 海葵目 12S 引物筛选','坐标为 Metridium senile NC_000933.1 的12S基因内位置。本次独立12S目标面板含'+str(len(panel))+'个具名物种代表。背景为'+str(len(back))+'条其他刺胞动物记录，各基因库独立取样。覆盖要求每端最多2个错配且各自末端3 nt兼容。预测目标产物100–500 bp，内部至少50 bp。严格 WC/GT 各≤10条、宽松 WC3/GT2 各≤50条为筛选门槛。普通与D共用起止坐标，D每端≤8变体，简并限于3′半段。新窗口先过背景门槛，再按覆盖数及完全匹配数选配。普通方案借鉴 MiFish 的保守侧翼设计思路，在目标模板和共识序列之间独立选配；D方案借鉴 Anth-28S 的谱系变异补偿思路，按同源模板变异生成简并候选，再按原始序列的成对覆盖、完全匹配与背景门槛比较。简并位置和简并度属于本项目约束。所有数字为计算预测。','|名称|参考跨度|可评价|覆盖|完全匹配|严格WC/GT|宽松WC3/GT2|通过背景门槛|','|---|---:|---:|---:|---:|---:|---:|---|']
for q in out:lines.append('| '+' | '.join([q['name'],q['span'],q['evaluable'],q['covered'],q['exact'],q['WC']+'/'+q['GT'],q['WC3']+'/'+q['GT2'],q['pass_gate']])+' |')
(p/'README.md').write_text('\n\n'.join(lines[:2])+'\n\n'+'\n'.join(lines[2:])+'\n');print('PASS',checks)
