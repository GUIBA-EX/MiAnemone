"""Evaluate target and background under identical length and mismatch models."""
from pathlib import Path
import csv
from collections import defaultdict
from independent import sites,comp
from verify import rows,fasta
ROOT=Path(__file__).resolve().parents[1]
def write(name,rr):
 with (ROOT/'revision'/name).open('w') as f:
  w=csv.DictWriter(f,rr[0],lineterminator='\n');w.writeheader();w.writerows(rr)
summary=[];lengths=[];raw=[]
for gene in ['12S','16S']:
 ps=rows(ROOT/gene/'candidate_results.csv');detail=rows(ROOT/gene/'species_results.csv')
 eligible={p['name']:{r['accession'] for r in detail if r['name']==p['name'] and r['evaluable']=='1'} for p in ps}
 for p in ps:
  for role in ['target','background']:
   lib=fasta(ROOT/gene/(role+'.fasta'));lib={k:s for k,s in lib.items() if role=='background' or k in eligible[p['name']]}
   for gt in [False,True]:
    counts=defaultdict(int)
    for acc,s in lib.items():
     pairs=[]
     for strand in [s,s.translate(comp)[::-1]]:
      fs=list(sites(strand,p['F'],False,gt,4));rs=list(sites(strand,p['R'],True,gt,4))
      for a,m,_ in fs:
       for b,n,_ in rs:
        span=b+len(p['R'])-a;gap=b-a-len(p['F'])
        if not 100<=span<=500 or gap<50:continue
        # Recalculate mismatch locations in primer orientation, GT as implemented in scanner.
        from independent import pattern
        ff,_=pattern(p['F'],False,gt,4);rr,_=pattern(p['R'],True,gt,4)
        df=[len(ff)-j-1 for j,c in enumerate(strand[a:a+len(ff)]) if c not in ff[j]]
        dr=[j for j,c in enumerate(strand[b:b+len(rr)]) if c not in rr[j]]
        pairs.append((df,dr,span))
     for k in range(5):
      for tail in [0,3,5]:
       hit=any(max(len(f),len(r))<=k and all(v>=tail for v in f+r) for f,r,span in pairs)
       counts[k,tail]+=hit
       if k in [2,3] and tail==3:
        raw.append(dict(gene=gene,pair=p['name'],role=role,model='GT' if gt else 'WC',accession=acc,max_mismatches=k,terminal=tail,hit=int(hit)))
     if role=='target' and not gt:
      ls=sorted({span for f,r,span in pairs if max(len(f),len(r))<=2 and all(v>=3 for v in f+r)})
      for span in ls:lengths.append(dict(gene=gene,pair=p['name'],accession=acc,product_bp=span,overlap_2x150=max(0,300-span),overlap_2x250=max(0,500-span),overlap_2x300=max(0,600-span)))
    for (k,tail),v in counts.items():summary.append(dict(gene=gene,pair=p['name'],role=role,model='GT' if gt else 'WC',max_mismatches=k,terminal_zero=tail,n=len(lib),matched=v))
  print(gene,p['name'],flush=True)
write('unified_models.csv',summary);write('unified_record_hits.csv',raw);write('product_lengths.csv',lengths)
