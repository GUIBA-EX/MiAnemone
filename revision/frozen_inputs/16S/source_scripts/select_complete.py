"""One best available 16S representative per named GenBank species."""
import csv,subprocess,tempfile,collections
from pathlib import Path
from Bio import SeqIO
p=Path(__file__).resolve().parents[1]
def read(n):return list(csv.DictReader((p/n).open()))
def write(n,rows):
 with (p/n).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
audit=read('accession_audit.csv');seq={r.id:str(r.seq) for r in SeqIO.parse(p/'all_eligible.fasta','fasta')};anchor='NC_000933.1'
rr=[x for x in audit if x['status']=='eligible' and x['named_binomial']=='1']
with tempfile.TemporaryDirectory() as tmp:
 t=Path(tmp);(t/'in.fa').write_text(''.join('>'+a+'\n'+seq[a]+'\n' for a in [anchor]+sorted(x['accession'] for x in rr if x['accession']!=anchor)))
 with (t/'al.fa').open('w') as out:subprocess.run(['mafft','--quiet','--adjustdirection','--retree','2','--maxiterate','0','--thread','4',str(t/'in.fa')],stdout=out,check=True)
 al={r.id.removeprefix('_R_'):str(r.seq).upper() for r in SeqIO.parse(t/'al.fa','fasta')}
cols=[i for i,b in enumerate(al[anchor]) if b in 'ACGT'];groups=collections.defaultdict(list);scores={}
for x in rr:
 a=x['accession']; called=sum(al[a][i] in 'ACGT' for i in cols);full=int('complete genome' in x['description'].lower() and x['method'].startswith('feature') and not set('<>')&set(x['location']) and called>=.9*len(cols))
 scores[a]=(full,called,int(x['length'])-int(x['ambiguous']),-int(x['ambiguous']),x['method'].startswith('feature'),a.startswith('NC_'),a)
 groups[x['organism']].append(x)
winners=[];ranks=[]
for name,rows in sorted(groups.items()):
 ranked=sorted(rows,key=lambda x:scores[x['accession']],reverse=True)
 for rank,x in enumerate(ranked,1):
  a=x['accession'];s=scores[a];x['representative']=int(rank==1)
  ranks.append(dict(organism=name,accession=a,rank=rank,selected=rank==1,records_for_species=len(rows),genome_supported_complete=s[0],called_reference_columns=s[1],reference_columns=len(cols),length=x['length'],ambiguous=x['ambiguous']))
 winners.append(ranked[0])
write('selection_ranking.csv',ranks);write('species_panel.csv',winners);write('accession_audit.csv',audit)
(p/'representatives.fasta').write_text(''.join(f">{x['accession']}\n{seq[x['accession']]}\n" for x in winners))
with (p/'alignment.fasta').open('w') as out:subprocess.run(['mafft','--quiet','--adjustdirection','--retree','2','--maxiterate','0','--thread','4',str(p/'representatives.fasta')],stdout=out,check=True)
