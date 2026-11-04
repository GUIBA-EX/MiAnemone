from pathlib import Path
import json,csv,collections,itertools,subprocess,shutil
from Bio import SeqIO
from chemistry import expand,thermal,deg,terminal
p=Path(__file__).resolve().parents[1];d=p/'12S';al={r.id.removeprefix('_R_'):str(r.seq).upper() for r in SeqIO.parse(d/'aligned.fasta','fasta')};cols=[i for i,b in enumerate(al['NC_000933.1']) if b!='-'];panel=list(csv.DictReader((d/'target_panel.csv').open()))
def write(path,rows):
 with path.open('w') as h:
  w=csv.DictWriter(h,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def score(where,qs):
 data={};eligible={}
 for q in qs:
  a,b,c,e=[cols[q[k]-1] for k in ['F_start','F_end','R_start','R_end']];ids={x['accession'] for x in panel if (s:=al[x['accession']]) and s[:a+1].strip('-') and s[e:].strip('-') and set(s[a:b+1]+s[c:e+1])<=set('ACGT-')};eligible[q['id']]=ids;data[q['id']]=dict(evaluable=len(ids),covered=0,exact=0,WC=0,GT=0,WC3=0,GT2=0)
 with (where/'independent_hits.tsv').open() as h:
  for x in csv.DictReader(h,delimiter='\t'):
   z=data[x['name']]
   if x['role']=='target':
    if x['accession'] in eligible[x['name']]:z['covered']+=int(x['le2last3']);z['exact']+=int(x['exact'])
   else:
    z[x['model']]+=int(x['le2last3']);z['WC3' if x['model']=='WC' else 'GT2']+=int(x['le3' if x['model']=='WC' else 'le2'])
 return [dict({k:v for k,v in q.items() if k!='thermals'},**data[q['id']]) for q in qs]
def passgate(q):return max(q['WC'],q['GT'])<=10 and max(q['WC3'],q['GT2'])<=50
