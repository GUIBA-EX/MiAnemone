"""Genus-blocked, fixed-window ordinary/D redesign. Endpoints are conditioned on full panel.
Candidate creation and thermodynamic gates use training genera only; background is not used to tune these diagnostic candidates.
"""
import hashlib,itertools,csv,collections
from pathlib import Path
from Bio.Seq import Seq
from verify import fasta,rows
from independent import flags
from review_chemistry import I,deg,terminal,thermal
ROOT=Path(__file__).resolve().parents[1];rev={frozenset(v):k for k,v in I.items() if k!='X'};out=[];choices=[]
def coverage(q,t):
 bad=[[i for i,(x,y) in enumerate(zip(a,b)) if y not in I[x]] for a,b in zip(q,t)]
 return all(len(z)<=2 and all(j<len(q[i])-3 for j in z) for i,z in enumerate(bad)),not any(bad)
for gene in ['12S','16S']:
 d=ROOT/gene;meta={r['accession']:r for r in rows(d/'target_panel.csv')};al={k.removeprefix('_R_'):v for k,v in fasta(d/'aligned.fasta').items()};co=[i for i,c in enumerate(al['NC_000933.1']) if c!='-'];raw=fasta(d/'target.fasta');ps=rows(d/'candidate_results.csv')
 for p in ps:
  if p['name'].endswith('D'):continue
  samples={}
  for k in meta:
   pair=[]
   for side in 'FR':
    t=al[k][co[int(p[side+'_start'])-1]:co[int(p[side+'_end'])-1]+1].replace('-','');pair.append(str(Seq(t).reverse_complement()) if side=='R' else t)
   if all(len(v)==len(p[side]) for side,v in zip('FR',pair)) and set(''.join(pair))<=set('ACGT'):samples[k]=tuple(pair)
  folds={k:int(hashlib.sha256(meta[k]['organism'].split()[0].encode()).hexdigest(),16)%5 for k in samples}
  for fold in range(5):
   train=[v for k,v in samples.items() if folds[k]!=fold];test=[k for k in samples if folds[k]==fold]
   if not train or not test:continue
   pools=[]
   for end in range(2):
    freq=collections.Counter(v[end] for v in train);L=len(next(iter(freq)));cons=''.join(max('ACGT',key=lambda b:sum(n for t,n in freq.items() if t[j]==b)) for j in range(L));ordinary={s for s in [*freq,cons] if terminal(s)}
    def rank(q):
     cov=ex=0
     for t,n in freq.items():
      bad=[i for i,(a,b) in enumerate(zip(q,t)) if b not in I[a]];cov+=n*(len(bad)<=2 and all(i<L-3 for i in bad));ex+=n*(not bad)
     return cov,ex,-deg(q),q
    mixes=set()
    for seed in sorted(ordinary,key=rank,reverse=True)[:12]:
     z=seed
     for _ in range(3):
      opts=[]
      for t in freq:
       if t[:L//2]!=seed[:L//2]:continue
       m=''.join(rev[frozenset(I[a]+b)] for a,b in zip(z,t))
       if m!=z and deg(m)<=8 and terminal(m):opts.append(m)
      if not opts:break
      mixes.update(opts);z=max(opts,key=rank)
    pools.append(sorted(ordinary,key=rank,reverse=True)[:6]+sorted(mixes,key=rank,reverse=True)[:6])
   candidates={False:[],True:[]}
   for q in itertools.product(*pools):
    ok,rep=thermal(*q)
    if not ok or min(min(r['F_tm_min'],r['R_tm_min']) for r in rep)<45:continue
    scores=[coverage(q,t) for t in train];rank=(sum(s[0] for s in scores),sum(s[1] for s in scores),-deg(q[0])*deg(q[1]),q);candidates[max(map(deg,q))>1].append(rank)
   for D,cc in candidates.items():
    q=max(cc)[-1] if cc else None
    choices.append(dict(gene=gene,window=p['name'],degenerate=D,fold=fold,training_n=len(train),test_n=len(test),eligible_candidates=len(cc),F=q[0] if q else '',R=q[1] if q else ''))
    for k in test:
     exact,covered,_=flags(raw[k],*q,True,False) if q else (False,False,False)
     out.append(dict(gene=gene,window=p['name'],degenerate=D,fold=fold,accession=k,genus=meta[k]['organism'].split()[0],candidate_available=q is not None,covered=int(covered),exact=int(exact)))
  print(gene,p['name'],flush=True)
for name,rr in [('redesign_holdout',out),('redesign_candidates',choices)]:
 with (ROOT/'revision'/(name+'.csv')).open('w') as h:
  w=csv.DictWriter(h,rr[0]);w.writeheader();w.writerows(rr)
