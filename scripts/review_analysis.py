"""Supplementary homologous-site diagnostics; paired raw-scan results remain separately identified."""
import csv,math,hashlib
from collections import Counter,defaultdict
from pathlib import Path
from Bio.Seq import Seq
from Bio.Data.IUPACData import ambiguous_dna_values as I
from verify import fasta,rows

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'revision';OUT.mkdir(exist_ok=True)
def write(name,data):
 if not data:return
 with (OUT/(name+'.csv')).open('w') as f:
  w=csv.DictWriter(f,data[0],lineterminator='\n');w.writeheader();w.writerows(data)
def wilson(k,n):
 if not n:return 0,0
 z=1.96;p=k/n;den=1+z*z/n;v=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))
 return (p+z*z/(2*n)-v)/den,(p+z*z/(2*n)+v)/den
out=defaultdict(list)
for gene in ['12S','16S']:
 d=ROOT/gene;al={k.removeprefix('_R_'):v for k,v in fasta(d/'aligned.fasta').items()};co=[i for i,c in enumerate(al['NC_000933.1']) if c!='-']
 ps=rows(d/'candidate_results.csv');meta={r['accession']:r for r in rows(d/'target_panel.csv')};det=rows(d/'species_results.csv');by={p['name']:[r for r in det if r['name']==p['name']] for p in ps}
 eligible={n:{r['accession'] for r in rr if r['evaluable']=='1'} for n,rr in by.items()};common=set.intersection(*eligible.values())
 resolution_common=set(common)
 for q in ps:
  left,right=co[int(q['F_end'])-1],co[int(q['R_start'])-1]
  resolution_common&={k for k in common if al[k][left+1:right].replace('-','') and set(al[k][left+1:right].replace('-',''))<=set('ACGT')}
 for j,c in enumerate(co):
  cs=Counter(al[k][c] for k in meta if al[k][c] in 'ACGT');n=sum(cs.values());entropy=-sum(v/n*math.log2(v/n) for v in cs.values()) if n else 0
  out['entropy'].append(dict(gene=gene,coordinate=j+1,n=n,entropy=entropy))
 for p in ps:
  name=p['name'];family=defaultdict(list)
  for r in by[name]:family[r['family']].append(r)
  for f,rr in family.items():
   out['family_coverage'].append(dict(gene=gene,pair=name,family=f,total=len(rr),evaluable=sum(int(r['evaluable']) for r in rr),covered=sum(int(r['covered']) for r in rr),exact=sum(int(r['exact']) for r in rr)))
  for label,ids in [('own',eligible[name]),('common',common)]:
   rr=[r for r in by[name] if r['accession'] in ids];n=len(rr);cov=sum(int(r['covered']) for r in rr);ex=sum(int(r['exact']) for r in rr);lo,hi=wilson(cov,n)
   out['panel_comparison'].append(dict(gene=gene,pair=name,panel=label,n=n,covered=cov,exact=ex,wilson_lower=lo,wilson_upper=hi))
  if name.endswith('D'):
   ordinary=next(r for r in ps if r['name']==name[:-1]);changes=[]
   for side in 'FR':
    for j,(a,b) in enumerate(zip(ordinary[side],p[side])):
     if not set(I[a])<=set(I[b]):changes.append(f'{side}{j+1}:{a}>{b}')
   dd={r['accession']:r for r in by[name]};oo={r['accession']:r for r in by[name[:-1]]}
   out['nesting'].append(dict(gene=gene,pair=name,nested=not changes,non_nested=';'.join(changes),ordinary_only=sum(int(oo[k]['covered']) and not int(dd[k]['covered']) for k in eligible[name]),D_only=sum(int(dd[k]['covered']) and not int(oo[k]['covered']) for k in eligible[name])))
  words={}; a,b,c,e=[co[int(p[k])-1] for k in ['F_start','F_end','R_start','R_end']]
  for k in eligible[name]:
   F=''.join(al[k][i] for i in co[int(p['F_start'])-1:int(p['F_end'])]);R=str(Seq(''.join(al[k][i] for i in co[int(p['R_start'])-1:int(p['R_end'])])).reverse_complement())
   clean=set(F+R)<=set('ACGT')
   words[k]=(F,R)
   for side,t in [('F',F),('R',R)]:
    for j,(q,z) in enumerate(zip(p[side],t)):
     out['mismatch_positions'].append(dict(gene=gene,pair=name,accession=k,species=meta[k]['organism'],side=side,distance_from_3prime=len(t)-j-1,primer=q,template=z,mismatch=int(z not in I[q]),gap=z=='-'))
   out['exclusions'].append(dict(gene=gene,pair=name,accession=k,raw_evaluable=True,complete_projected_sites=clean))
  for maxmm in range(4):
   for tail in [3,5]:
    cov=0;penalties=[]
    for k,(F,R) in words.items():
     if set(F+R)-set('ACGT'):continue
     bad=[[len(t)-j-1 for j,(q,z) in enumerate(zip(p[s],t)) if z not in I[q]] for s,t in [('F',F),('R',R)]]
     cov+=all(len(v)<=maxmm and all(x>=tail for x in v) for v in bad)
     penalties.append(sum(sum(1+4*math.exp(-x/3) for x in v) for v in bad))
    out['site_sensitivity'].append(dict(gene=gene,pair=name,mismatches=maxmm,terminal_zero=tail,n=len(penalties),covered=cov,mean_positional_penalty=sum(penalties)/len(penalties) if penalties else 0))
  # Reference insert uniqueness conditioned on own and common evaluability.
  for panel,ids in [('own',eligible[name]),('common',resolution_common)]:
   ins={k:al[k][b+1:c].replace('-','') for k in ids};ins={k:v for k,v in ins.items() if v and set(v)<=set('ACGT')};counts=Counter(ins.values());h={r['accession']:int(r['covered']) for r in by[name]}
   out['joint_resolution'].append(dict(gene=gene,pair=name,panel=panel,n=len(ins),unique=sum(counts[v]==1 for v in ins.values()),covered_unique=sum(counts[v]==1 and h[k] for k,v in ins.items())))
 # Five deterministic genus folds, reselect the observed paired motif by training exact then compatible.
 for p in ps:
  if p['name'].endswith('D'):continue
  name=p['name'];samples={}
  for k in eligible[name]:
   F=al[k][co[int(p['F_start'])-1]:co[int(p['F_end'])-1]+1].replace('-','');R=str(Seq(al[k][co[int(p['R_start'])-1]:co[int(p['R_end'])-1]+1].replace('-','')).reverse_complement())
   if len(F)==len(p['F']) and len(R)==len(p['R']) and set(F+R)<=set('ACGT'):samples[k]=(F,R)
  folds={k:int(hashlib.sha256(meta[k]['organism'].split()[0].encode()).hexdigest(),16)%5 for k in samples}
  for fold in range(5):
   train=[v for k,v in samples.items() if folds[k]!=fold];test=[k for k in samples if folds[k]==fold]
   if not train or not test:continue
   candidates=Counter(train)
   def score(q):
    exact=candidates[q];cov=sum(all(sum(x!=y for x,y in zip(a,b))<=2 and a[-3:]==b[-3:] for a,b in zip(q,t)) for t in train);return exact,cov,q
   q=max(candidates,key=score)
   for k in test:
    v=samples[k];out['genus_holdout'].append(dict(gene=gene,window=name,fold=fold,accession=k,genus=meta[k]['organism'].split()[0],F=q[0],R=q[1],exact=q==v,covered=all(sum(x!=y for x,y in zip(a,b))<=2 and a[-3:]==b[-3:] for a,b in zip(q,v))))
for name,data in out.items():write(name,data)
print({k:len(v) for k,v in out.items()})
