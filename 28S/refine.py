from pathlib import Path
import csv,json,collections,itertools
from Bio import SeqIO
from Bio.Seq import Seq
from chemistry import I,deg,terminal,thermal,expand
r=Path(__file__).resolve().parent;old={x['name']:x for x in csv.DictReader((r/'prior_candidates.csv').open())};al={x.id.removeprefix('_R_'):str(x.seq).upper() for x in SeqIO.parse(r/'target.aligned.fasta','fasta')};ref=json.loads((r/'reference.json').read_text())['accession'];co=[i for i,b in enumerate(al[ref]) if b!='-'];inv={frozenset(v):k for k,v in I.items() if k!='X'};audit=[];out=[];chem=[]
def masks(q,ts):
 cv=ex=0
 for i,t in enumerate(ts):
  if len(t)!=len(q) or set(t)-set('ACGT'):continue
  bad=[j for j,(a,b) in enumerate(zip(q,t)) if b not in I[a]]
  if not bad:ex|=1<<i
  if len(bad)<=2 and all(j<len(q)-3 for j in bad):cv|=1<<i
 return cv,ex
for name,source in [('XS','S1'),('S','S3'),('M','M2'),('L','L1')]:
 p=old[source];d=old[source+'D'];pools={};ordinary={}
 for side in 'FR':
  a,b=[co[int(p[side+'_'+z])-1] for z in ['start','end']];ts=[s[a:b+1].replace('-','') for s in al.values()]
  if side=='R':ts=[str(Seq(t).reverse_complement()) for t in ts]
  good=[t for t in ts if len(t)==len(p[side]) and set(t)<=set('ACGT')];ct=[collections.Counter(v) for v in zip(*good)];major=''.join(c.most_common(1)[0][0] for c in ct);cons=''
  for j,c in enumerate(ct):
   chosen=[];n=0
   for base,num in c.most_common():
    chosen.append(base);n+=num
    if n>=.95*len(good):break
   code=inv[frozenset(chosen)];cons+=code;audit.append(dict(window=name,side=side,position=j+1,n=len(good),A=c['A'],C=c['C'],G=c['G'],T=c['T'],majority=major[j],consensus95=code))
  seqs={p[side],d[side],major};front={major}
  for step in range(3):
   nxt=set()
   for q in front:
    for j,c in enumerate(cons):
     for b in I[c]:
      if b not in I[q[j]]:
       z=q[:j]+inv[frozenset(I[q[j]]+b)]+q[j+1:]
       if deg(z)<=8 and terminal(z):nxt.add(z)
   front=set(sorted(nxt,key=lambda q:tuple(x.bit_count() for x in masks(q,ts)),reverse=True)[:12]);seqs|=front
  if deg(cons)<=8:seqs.add(cons)
  pools[side]=[(q,*masks(q,ts)) for q in seqs if terminal(q)]
  concretes={p[side],major}|set(t for t,n in collections.Counter(good).most_common(12))|{s for q in seqs for s in expand(q) if deg(q)<=8}
  ordinary[side]=[(q,*masks(q,ts)) for q in concretes if terminal(q)]
 for mode,opts in [('ordinary',ordinary),('D',pools)]:
  pairs=sorted((((fc&rc).bit_count(),(fe&re).bit_count(),F,R) for F,fc,fe in opts['F'] for R,rc,re in opts['R'] if mode=='ordinary' or deg(F)*deg(R)>1),reverse=True);valid=[]
  for cv,ex,F,R in pairs:
   ok,rr=thermal(F,R);ok=ok and all(min(z['F_tm_min'],z['R_tm_min'])>=45 for z in rr)
   if ok:valid.append((cv,ex,F,R))
   if len(valid)>=30:break
  if not valid:raise ValueError((name,mode))
  best=max(v[0] for v in valid);valid=sorted([v for v in valid if v[0]>=best-2],key=lambda v:(v[1],v[0]),reverse=True)
  for j,(cv,ex,F,R) in enumerate(valid[:4]):
   z={k:p[k] for k in ['F_start','F_end','R_start','R_end','span']};z.update(name=f'{name}trial{j+1}'+('D' if mode=='D' else ''),F=F,R=R,deg_F=deg(F),deg_R=deg(R),thermal_pass=True);out.append(z)
 # preserve exact prior alternatives
 for suffix in ['', 'D']:
  z=dict(old[source+suffix]);z['name']=name+'prior'+suffix;out.append(z)
 print(name,'prepared',flush=True)
for fn,rr in [('consensus95_audit.csv',audit),('candidates.csv',out)]:
 with (r/fn).open('w') as h:w=csv.DictWriter(h,rr[0]);w.writeheader();w.writerows(rr)
(r/'independent_primers.tsv').write_text(''.join(f"{x['name']}\t{x['F']}\t{x['R']}\n" for x in out))
