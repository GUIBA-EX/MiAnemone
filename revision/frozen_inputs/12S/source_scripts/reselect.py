"""Rescreen fixed homologous windows with independent concrete and IUPAC pools."""
from choose import p,d,al,cols,panel,score,passgate,write
from chemistry import I,deg,terminal,thermal
from Bio.Seq import Seq
import shutil,csv,collections,itertools,subprocess,tempfile,json
draft=score(d,json.loads((d/'candidates.json').read_text()));write(p/'screened.csv',draft);previous=[]
for group in ('XS','S','M','L','XL'):
 pool=[q for q in draft if q['group']==group and passgate(q)]
 if pool:
  q=max(pool,key=lambda q:(q['covered'],q['exact'],-q['GT']));previous.append(dict(q,id=group+'D'))

rev={frozenset(v):k for k,v in I.items() if k!='X'}
probes=[];supports=[]
for group in [q['group'] for q in previous]:
 ref=next(q for q in previous if q['id']==group+'D')
 base={k:int(ref[k]) for k in ['F_start','F_end','R_start','R_end','span']}
 if group=='fixed':probes.append(dict(base,group=group,id='H',F='AAGWCCCCATGRARCTTTACTG',R='CCCCAGRGTAACTTTTATYTGY',deg_F=8,deg_R=8))
 sites={}
 for side in 'FR':
  a,b=base[side+'_start'],base[side+'_end'];L=b-a+1;freq=collections.Counter()
  for row in panel:
   t=al[row['accession']][cols[a-1]:cols[b-1]+1].replace('-','')
   if len(t)==L and set(t)<=set('ACGT'):freq[str(Seq(t).reverse_complement()) if side=='R' else t]+=1
  consensus=''.join(max('ACGT',key=lambda c:sum(v for t,v in freq.items() if t[i]==c)) for i in range(L));ordinary={t for t in [*freq,consensus] if terminal(t)}
  def rank(q):
   compatible=exact=0
   for t,n in freq.items():
    bad=[i for i,(x,y) in enumerate(zip(q,t)) if y not in I[x]]
    compatible+=n*(len(bad)<=2 and all(i<L-3 for i in bad));exact+=n*(not bad)
   return compatible,exact,-deg(q),q
  mix=set()
  for seed in sorted(ordinary,key=rank,reverse=True)[:12]:
   z=seed
   for _ in range(3):
    opts=[]
    for t in freq:
     if t[:L//2]!=seed[:L//2]:continue
     m=''.join(rev[frozenset(I[x]+y)] for x,y in zip(z,t))
     if m!=z and deg(m)<=8 and terminal(m):opts.append(m)
    if not opts:break
    mix.update(opts);z=max(opts,key=rank)
  old=ref[side]
  if terminal(old) and all(set(I[x]) <= {t[i] for t in freq} for i,x in enumerate(old)):mix.add(old)
  top=lambda v:sorted(v,key=rank,reverse=True)[:6]
  sites[side]=top(ordinary)+top({t for t in mix if deg(t)>1})
  for t,n in sorted(freq.items()):supports.append(dict(window=group,side=side,template=t,species=n))
 for i,(f,r) in enumerate(itertools.product(sites['F'],sites['R'])):
  ok,details=thermal(f,r)
  if ok and min(min(t['F_tm_min'],t['R_tm_min']) for t in details)>=45:probes.append(dict(base,group=group,id=f'{group}_{i}',F=f,R=r,deg_F=deg(f),deg_R=deg(r)))
write(p/'template_support.csv',supports)
v=p/'verify';v.mkdir(exist_ok=True);shutil.copy(d/'target.fasta',v/'target.fasta');shutil.copy(d/'background.fasta',v/'background.fasta');(v/'independent_primers.tsv').write_text(''.join('{id}\t{F}\t{R}\n'.format(**q) for q in probes))
with tempfile.TemporaryDirectory() as tmp,(v/'scan.log').open('w') as h:
 exe=tmp+'/scan';subprocess.run(['c++','-O3','-std=c++17',str(p/'scripts/raw_scan.cpp'),'-o',exe],check=True);subprocess.run([exe,str(v)],stderr=h,check=True)
rs=score(v,probes);write(p/'method_screen.csv',rs);res=[]
for group in [q['group'] for q in previous]:
 for D in (False,True):
  eligible=[q for q in rs if q['group']==group and (max(q['deg_F'],q['deg_R'])>1)==D]
  if group!='fixed':eligible=[q for q in eligible if passgate(q)]
  if not eligible:raise ValueError((group,D))
  q=max(eligible,key=lambda q:(q['covered'],q['exact'],-q['GT'],-q['WC3']))
  res.append(q)
for q in res:q['pass_gate']=passgate(q)
write(p/'comparison.csv',res)
