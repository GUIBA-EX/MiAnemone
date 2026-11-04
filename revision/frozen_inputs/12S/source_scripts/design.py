"""Screen fixed reference windows independently for one gene. Run after MAFFT."""
from pathlib import Path
import sys,csv,json,collections,itertools
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from chemistry import I,deg,terminal,thermal,oligo_thermal
root=Path(__file__).resolve().parents[1];d=root/sys.argv[1];anchor=(d/'anchor.txt').read_text().strip();al={r.id.removeprefix('_R_'):str(r.seq).upper() for r in SeqIO.parse(d/'aligned.fasta','fasta')};cols=[i for i,c in enumerate(al[anchor]) if c!='-'];panel=list(csv.DictReader((d/'target_panel.csv').open()));ids=[r['accession'] for r in panel];N=len(ids);bits=lambda a: int.from_bytes(np.packbits(a,bitorder='little').tobytes(),'little')
rev={frozenset(v):k for k,v in I.items() if k!='X'};masks={q:sum(1<<'ACGT'.index(b) for b in v) for q,v in I.items() if set(v)<=set('ACGT')}
sites={'F':[],'R':[]}
for side in sites:
 for a in range(0,len(cols)-22,3):
  local=[]
  for L in [20,22]:
   ts=[al[k][cols[a]:cols[a+L-1]+1].replace('-','') for k in ids];ts=[(str(Seq(t).reverse_complement()) if side=='R' else t) if len(t)==L and set(t)<=set('ACGT') else '-'*L for t in ts]
   freq=collections.Counter(t for t in ts if '-' not in t)
   if sum(freq.values())<max(10,N*.10):continue
   mat=np.array([[masks.get(c,0) for c in t] for t in ts],dtype='uint8');present=(mat!=0).all(axis=1)
   consensus=''.join(max('ACGT',key=lambda b:sum(n for t,n in freq.items() if t[j]==b)) for j in range(L));cand=set()
   for start in [consensus]+[t for t,n in freq.most_common(5)]:
    cand.add(start);q=start
    choices=sorted([(sum(n for t,n in freq.items() if t[j]==b),j,b) for j in range(L//2,L) for b in 'ACGT' if b!=start[j]],reverse=True)
    for n,j,b in choices:
     if n<.04*sum(freq.values()):break
     z=q[:j]+rev[frozenset(I[q[j]]+b)]+q[j+1:]
     if deg(z)<=8 and terminal(z):q=z;cand.add(q)
   for q in sorted(cand):
    if not terminal(q):continue
    bad=(mat&np.array([masks[c] for c in q]))==0;cov=bits(present & (bad.sum(axis=1)<=2)&~bad[:,-3:].any(axis=1));exact=bits(present & ~bad.any(axis=1))
    if cov:local.append((a,L,q,cov,exact))
  if local:
   score=lambda x:(x[3].bit_count(),x[4].bit_count(),-deg(x[2]))
   if False:sites[side]+=sorted([x for x in local if x[1]==22],key=score,reverse=True)[:35]
   sites[side]+=sorted(local,key=score,reverse=True)[:3]+sorted([x for x in local if deg(x[2])==1],key=score,reverse=True)[:1]
pairs=[]
for f,r in itertools.product(sites['F'],sites['R']):
 span=r[0]+r[1]-f[0]
 if 150<=span<=470:
  sc=(f[3]&r[3]).bit_count();ex=(f[4]&r[4]).bit_count()
  if sc:pairs.append((sc,ex,f,r,span))
# Keep diverse coordinate bins before raw background screening.
bybin=collections.defaultdict(list)
for row in sorted(pairs,key=lambda z:(z[0],z[1]),reverse=True):
 sc,ex,f,r,span=row;k=(f[0]//25,r[0]//25)
 if len(bybin[k])<3:bybin[k].append(row)
short=sorted([x for rr in bybin.values() for x in rr],key=lambda z:(z[0],z[1]),reverse=True)
if d.name=='COI':
 # Screen high-coverage windows plus coordinate-diverse windows across the gene.
 strata=sorted(bybin)
 order=[bybin[k][0] for k in strata[::max(1,len(strata)//450)]]
 short=short[:160]+order+short
fixed=[x for x in pairs if x[2][:2]==(1496,22) and x[3][:2]==(1808,22)];short=sorted(fixed,key=lambda z:(z[0],z[1]),reverse=True)+short
def band(f,r,span):
 
 return next(b for b,n in [('XS',185),('S',210),('M',250),('L',350),('XL',470)] if span<=n)
chosen=[];seen=set();counts=collections.Counter()
for sc,ex,f,r,span in short:
 k=f[2],r[2]
 group=band(f,r,span)
 if max(deg(f[2]),deg(r[2]))==1 or counts[group]>=90 or k in seen:continue
 seen.add(k)
 at,bt=oligo_thermal(k[0])[0],oligo_thermal(k[1])[0]
 if min(at[0],bt[0])<45 or max(at[2],bt[2])>=55:continue
 ok,th=thermal(*k)
 if not ok:continue
 counts[group]+=1
 chosen.append(dict(group=group,id='P'+str(len(chosen)+1),F=f[2],R=r[2],F_start=f[0]+1,F_end=f[0]+f[1],R_start=r[0]+1,R_end=r[0]+r[1],span=span,deg_F=deg(f[2]),deg_R=deg(r[2]),thermals=th))

if not chosen:raise RuntimeError('No thermally eligible candidate')
(d/'candidates.json').write_text(json.dumps(chosen,indent=2));(d/'independent_primers.tsv').write_text(''.join(f"{p['id']}\t{p['F']}\t{p['R']}\n" for p in chosen))
print(d.name,len(pairs),len(chosen),'screened',flush=True)
