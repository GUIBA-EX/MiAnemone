import sys,csv,json,shutil
from pathlib import Path
r=Path(__file__).resolve().parents[1];sys.path.insert(0,str(r/'scripts'))
from verify import rows,fasta
from independent import flags
from review_chemistry import thermal

def write(p,rr):
 with p.open('w') as f:
  w=csv.DictWriter(f,rr[0],lineterminator='\n');w.writeheader();w.writerows(rr)
archive=r/'revision/prior_XSD';archive.mkdir(exist_ok=True)
for n in ['factorial_controls.csv','factorial_filters.csv']:
 if not (archive/n).exists():shutil.copy2(r/'revision'/n,archive/n)
ps=rows(r/'12S/candidate_results.csv');q=next(x for x in ps if x['name']=='XSD')
if not (archive/'candidate_results.csv').exists():write(archive/'candidate_results.csv',[q])
q.update(F='CGAAAGGTAAAATCTTCCGCCT',R='CGACGGCCATGCARTACCDG',id='XS_nested',deg_F='1',deg_R='6')
det=rows(r/'12S/species_results.csv');raw=fasta(r/'12S/target.fasta')
for x in det:
 if x['name']=='XSD':
  ex,cv,_=flags(raw[x['accession']],q['F'],q['R'],True,False);x['covered']=int(cv and int(x['evaluable']));x['exact']=int(ex and int(x['evaluable']))
for k in ['covered','exact']:q[k]=sum(int(x[k]) for x in det if x['name']=='XSD')
hits=[x for x in rows(r/'12S/background_hits.csv') if x['name']!='XSD'];meta={x['accession']:x for x in rows(r/'12S/background_panel.csv')}
for k in ['WC','GT','WC3','GT2']:q[k]=0
for acc,s in fasta(r/'12S/background.fasta').items():
 for model,broad in [('WC','WC3'),('GT','GT2')]:
  _,a,b=flags(s,q['F'],q['R'],False,model=='GT');q[model]+=int(a);q[broad]+=int(b)
  if a or b:hits.append(dict(name='XSD',accession=acc,organism=meta[acc]['organism'],lineage=meta[acc]['lineage'],model=model,strict=int(a),relaxed=int(b)))
write(r/'12S/candidate_results.csv',ps);write(r/'12S/species_results.csv',det);write(r/'12S/background_hits.csv',hits)
allp=rows(r/'all_results.csv')
for x in allp:
 if x['gene']=='12S' and x['name']=='XSD':x.update({k:v for k,v in q.items() if k in x})
write(r/'all_results.csv',allp)
for gene in ['12S','16S']:
 gg=rows(r/gene/'candidate_results.csv');(r/gene/'primers.fasta').write_text(''.join(f'>{gene}_{x["name"]}_{s}\n{x[s]}\n' for x in gg for s in 'FR'))
(r/'primers.fasta').write_text(''.join(f'>{x["gene"]}_{x["name"]}_{s}\n{x[s]}\n' for x in allp for s in 'FR'))
(r/'revision/XSD_replacement.json').write_text(json.dumps(dict(selected=q,thermal=thermal(q['F'],q['R'])),indent=2))
print(q)
