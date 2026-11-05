from pathlib import Path
import csv,json,math,hashlib
from Bio.Data.IUPACData import ambiguous_dna_values as I
r=Path(__file__).resolve().parent
read=lambda f:list(csv.DictReader((r/f).open()))
rows=read('comparison.csv');res={x['name']:x for x in read('rescue_summary.csv')};selected=[];changes=[];mapping={};details=[]
for n in ['XS','S','M','L']:
 for d in [False,True]:
  name=n+('D' if d else '');prior=next(x for x in rows if x['name']==n+'prior'+('D' if d else ''))
  pool=[x for x in rows if x['name'].startswith(n+'trial') and x['name'].endswith('D')==d]+[prior]
  pool=[x for x in pool if int(x['bg_WC'])<=int(prior['bg_WC']) and int(x['bg_GT'])<=int(prior['bg_GT']) and int(x['covered'])>=int(prior['covered'])-2]
  best=max(int(x['covered']) for x in pool);pool=[x for x in pool if int(x['covered'])>=best-2]
  chosen=max(pool,key=lambda x:(int(x['exact']),int(x['covered']),-int(x['deg_F'])-int(x['deg_R'])))
  z=dict(chosen);mapping[z['name']]=name;source=z.pop('name');z={'name':name,**z};selected.append(z)
  changes.append(dict(name=name,coverage_before=prior['covered'],coverage_after=z['covered'],exact_before=prior['exact'],exact_after=z['exact'],bg_WC_before=prior['bg_WC'],bg_WC_after=z['bg_WC'],bg_GT_before=prior['bg_GT'],bg_GT_after=z['bg_GT']))
  q=res[source];details.append(dict(name=name,core_pairs=q['unresolved_across_both_markers'],core_resolved_covered=q['core_resolved_and_covered'],any_window_resolved_covered=q['resolved_and_both_covered']))
for fn,rr in [('selected.csv',selected),('changes.csv',changes),('selected_rescue.csv',details)]:
 with (r/fn).open('w') as h:w=csv.DictWriter(h,rr[0]);w.writeheader();w.writerows(rr)
for f in ['species_results.csv','target_by_family.csv','background_by_group.csv']:
 rr=[dict(x,name=mapping[x['name']]) for x in read(f) if x['name'] in mapping]
 with (r/('selected_'+f)).open('w') as h:w=csv.DictWriter(h,rr[0]);w.writeheader();w.writerows(rr)
(r/'primers.fasta').write_text(''.join(f">28S_{x['name']}_{side}\n{x[side]}\n" for x in selected for side in 'FR'))
for i in range(0,8,2):
 a,b=selected[i:i+2];assert all(a[k]==b[k] for k in ['F_start','F_end','R_start','R_end'])
 for q in [a,b]:
  for side in 'FR':assert math.prod(len(I[c]) for c in q[side])<=8 and len(q[side])==int(q[side+'_end'])-int(q[side+'_start'])+1
v=json.loads((r/'verification.json').read_text());assert not v['disagreements']
lines=['# MiAnemone 28S：XS、S、M、L 复评','','## 方法','','基于221个海葵物种代表及448个其他六放珊瑚物种背景。参考 Anth-28S 先计算每个位点覆盖至少95%观测碱基的共识，再按真实变异生成低简并候选。95%是逐位点阈值，不是整对引物覆盖率；每条简并度≤8，故不保证完整实现全部95%共识。','', '普通引物比较共识、常见实际序列及D的具体变体；先通过末端和热力学门槛，再在近最大覆盖候选中优先完全匹配。独立重新扫描全部原始序列。保留旧方案作对照，最终要求严格WC/GT背景不增加、覆盖下降不超过2个物种，再按完全匹配、覆盖和较低简并度选择。此次为固定窗口再筛选，不是重新搜索全基因。','','所有具体变体在50 nM/变体及200 nM/方向下检查：结构Tm<55°C、最低结合Tm≥45°C、F/R Tm范围中点差≤5°C。普通和D的4个端点完全一致。','','[Anth-28S 原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC11606329/) · [95%共识逐碱基证据](consensus95_audit.csv)','','## 保留的方案','','|名称|参考跨度|F（5′–3′）|R（5′–3′）|F/R简并度|覆盖|完全匹配|严格WC/GT|宽松WC3/GT2|','|---|---:|---|---|---|---|---|---|---|']
for x in selected:lines.append(f"|{x['name']}|{x['span']}|{x['F']}|{x['R']}|{x['deg_F']}/{x['deg_R']}|{x['covered']}/{x['n']}|{x['exact']}/{x['n']}|{x['bg_WC']}/{x['bg_GT']}|{x['bg_WC3']}/{x['bg_GT2']}|")
lines+=['','## 补充分辨率','','以下112对物种在所有可比较的12S与16S窗口均共享，且两个基因均有可比较窗口。计数要求28S序列不同且两个物种均被覆盖；缺失28S不计作分辨失败。','','|方案|区分且两物种覆盖的物种对|','|---|---:|']
for x in details:lines.append(f"|{x['name']}|{x['core_resolved_covered']}|")
lines+=['','## 结论','','L及LD窗口主要提供分辨率增益，XS/XSD主要提供短片段与匹配率比较。严格模型零背景命中仍高度依赖末端约束，宽松模型有大量其他六放珊瑚命中；未证明实验特异性。普通方案较低的完全匹配率如未改善，保留为同坐标比较，不等同于推荐。核与线粒体参考按物种名关联，不代表相同个体。','','[变化表](changes.csv) · [序列FASTA](primers.fasta) · [各科覆盖](selected_target_by_family.csv) · [背景类群](selected_background_by_group.csv) · [独立复核](verification.json)']
(r/'report.md').write_text('\n'.join(lines)+'\n');(r/'validation.json').write_text(json.dumps(dict(status='PASS',names=[x['name'] for x in selected],checks=v['checks'],disagreements=0,same_coordinates=True,max_degeneracy=8),indent=2)+'\n');print('PASS',v['checks']);print('\n'.join(str(x) for x in changes))
