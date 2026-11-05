import urllib.request as u,urllib.parse as p,time,json,xml.etree.ElementTree as E
from pathlib import Path
out=Path(__file__).resolve().parent/'data'
def call(api,data):
 for attempt in range(6):
  try:
   time.sleep(.4);return u.urlopen('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'+api+'?'+p.urlencode(data),timeout=60).read()
  except Exception:
   if attempt==5:raise
   time.sleep(2**attempt)
for g in ['Actiniaria','Scleractinia','Corallimorpharia','Zoantharia','Antipatharia','Ceriantharia']:
 term=f'{g}[Organism] AND (28S[Title] OR "large subunit ribosomal RNA"[Title] OR "large subunit ribosomal DNA"[Title]) NOT mitochondrial[Title] NOT mitochondrion[Title] AND 150:15000[SLEN]'
 q=out/(g+'_search.json')
 if not q.exists():q.write_bytes(call('esearch.fcgi',dict(db='nuccore',term=term,retmax=100000,retmode='json')))
 a=json.loads(q.read_text())['esearchresult'];ids=a['idlist'];assert len(ids)==int(a['count']),(g,a['count'])
 (out/(g+'_query.txt')).write_text(term+'\nRetrieved 2026-11-05\n')
 for j in range(0,len(ids),100):
  f=out/f'{g}_{j:06d}.xml'
  if not f.exists():f.write_bytes(call('efetch.fcgi',dict(db='nuccore',id=','.join(ids[j:j+100]),rettype='gb',retmode='xml')))
  print(g,j,len(ids),flush=True)
