from pathlib import Path
import csv,json,sys
from Bio import SeqIO
r=Path(__file__).resolve().parent;sys.path.insert(0,str(r))
from independent import flags
read=lambda f:list(csv.DictReader(f.open()))
ps=read(r/'candidates.csv');hh=list(csv.DictReader((r/'independent_hits.tsv').open(),delimiter='\t'));h={(x['name'],x['role'],x['model'],x['accession']):x for x in hh};checks=0;bad=[]
for role in ['target','background']:
 for rec in SeqIO.parse(r/(role+'.fasta'),'fasta'):
  for p in ps:
   for gt in [False,True]:
    exact,covered,_=flags(str(rec.seq).upper(),p['F'],p['R'],True,gt);x=h.get((p['name'],role,'GT' if gt else 'WC',rec.id),{});checks+=2
    if (int(exact),int(covered))!=(int(x.get('exact',0)),int(x.get('le2last3',0))):bad.append([role,p['name'],rec.id,gt,exact,covered,x])
 print(role,'checked',checks,flush=True)
(r/'verification.json').write_text(json.dumps(dict(checks=checks,disagreements=bad,models='Strict WC and GT, same 100–500 bp product geometry'),indent=2)+'\n');assert not bad
