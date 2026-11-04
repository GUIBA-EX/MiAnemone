from pathlib import Path
import urllib.request as req,urllib.parse as enc,json,gzip,time
p=Path(__file__).resolve().parents[1]/'downloads'
def api(name,**args):
 for i in range(5):
  try:
   with req.urlopen('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'+name+'.fcgi?'+enc.urlencode(dict(db='nuccore',**args)),timeout=120) as f:return f.read()
  except Exception:
   if i==4:raise
   time.sleep(i+1)
q='Actiniaria[Organism] AND (16S[All Fields] OR rrnL[Gene Name] OR rnl[Gene Name] OR large subunit ribosomal RNA[Title] OR (mitochondrion[Filter] AND 10000:100000[SLEN]) OR (mitochondrial[Title] AND genome[Title]))'
ids=[]
while True:
 v=json.loads(api('esearch',term=q,retmode='json',retmax=5000,retstart=len(ids)))['esearchresult'];ids+=v['idlist']
 if len(ids)>=int(v['count']):break
 assert v['idlist'];time.sleep(.4)
assert len(ids)==len(set(ids))==int(v['count']);(p/'search.json').write_text(json.dumps(dict(query=q,translation=v.get('querytranslation'),count=len(ids),uids=ids),indent=2));print(len(ids),flush=True)
for i in range(0,len(ids),150):
 chunk=ids[i:i+150];raw=api('efetch',id=','.join(chunk),rettype='gb',retmode='text');assert raw.count(b'\n//')==len(chunk);(p/f'{i:06}.gb.gz').write_bytes(gzip.compress(raw));print(i+len(chunk),flush=True);time.sleep(.4)
