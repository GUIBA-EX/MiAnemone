import collections,csv,json,hashlib,zipfile
from choose import p,d,panel,al,cols,write
from chemistry import thermal,deg,I,terminal
from Bio import SeqIO
rr=list(csv.DictReader((p/'candidate_results.csv').open()));val=json.loads((p/'validation.json').read_text());qc=[];coords=[];details=[]
assert val['disagreements']==0 and len({r['name'] for r in rr})==len(rr)
fa={r.id:str(r.seq) for r in SeqIO.parse(p/'primers.fasta','fasta')}
for band in ['XS','S','M','L','XL']:
 pairs=[r for r in rr if r['group']==band]
 if not pairs:continue
 a,b=pairs;assert a['name']==band and b['name']==band+'D'
 keys=['F_start','F_end','R_start','R_end','span'];assert all(a[k]==b[k] for k in keys);coords.append(dict(window=band,reference='NC_000933.1',**{k:a[k] for k in keys}))
 for r in pairs:
  for side in 'FR':
   q=r[side];assert fa[r['name']+'_'+side]==q and deg(q)==int(r['deg_'+side])<=8 and terminal(q) and set(q[:len(q)//2])<=set('ACGT') and len(q)==int(r[side+'_end'])-int(r[side+'_start'])+1
  assert max(int(r[k]) for k in ['deg_F','deg_R'])>1 if r is b else max(int(r[k]) for k in ['deg_F','deg_R'])==1
  ok,t=thermal(r['F'],r['R']);assert ok and min(min(v['F_tm_min'],v['R_tm_min']) for v in t)>=45 and r['pass_gate']=='True'
  qc.append(dict(name=r['name'],deg_F=r['deg_F'],deg_R=r['deg_R'],max_structure_Tm=max(v['structure_max'] for v in t),status='PASS'))
 a0,b0,c,e=[cols[int(a[k])-1] for k in keys[:4]];bars=[]
 for x in panel:
  s=al[x['accession']];ins=s[b0+1:c].replace('-','')
  if s[:a0+1].strip('-') and s[e:].strip('-') and set(s[a0:b0+1]+s[c:e+1])<=set('ACGT-') and ins and set(ins)<=set('ACGT'):bars.append(ins)
 details.append(dict(window=band,unambiguous=len(bars),unique=sum(v==1 for v in collections.Counter(bars).values())))
write(p/'coordinates.csv',coords);write(p/'thermal_check.csv',qc);write(p/'barcode_identity.csv',details)
with (p/'README.md').open('a') as f:f.write('\n目标库来源为本轮NCBI检索，按注释和明确碱基数量为具名物种选代表；背景为已有刺胞动物GenBank库中可提取12S的具名物种，包含216种石珊瑚和8种珊瑚形目。各窗口分母按两端可评价性单独计算。\n\n[候选序列](primers.fasta) · [逐对指标](candidate_results.csv) · [参考结合坐标](coordinates.csv) · [内部条码唯一性](barcode_identity.csv) · [独立复核](validation.json)\n')
with zipfile.ZipFile(p.parent/'Actiniaria_12S.zip','w',zipfile.ZIP_DEFLATED) as z:
 for f in sorted(p.rglob('*')):
  if f.is_file() and '__pycache__' not in f.parts and f.suffix!='.log':z.write(f,str(f.relative_to(p.parent)))
