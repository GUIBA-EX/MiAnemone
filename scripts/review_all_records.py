"""Reconstruct named-species record panels from frozen archives and align all records."""
import csv,gzip,re,subprocess,json
from pathlib import Path
from Bio import SeqIO
from verify import rows,fasta
ROOT=Path(__file__).resolve().parents[1];WORK=ROOT/'revision/frozen_inputs'
for gene in ['12S','16S']:
 out=ROOT/'revision'/('all_records_'+gene);out.mkdir(exist_ok=True);seq={};meta=[]
 if gene=='16S':
  seq0=fasta(WORK/'16S/all_eligible.fasta')
  for r in rows(WORK/'16S/accession_audit.csv'):
   if r['status']=='eligible' and r['named_binomial']=='1':seq[r['accession']]=seq0[r['accession']];meta.append(r)
 else:
  audit={r['accession']:r for r in rows(WORK/'12S/audit.csv') if r['role']=='target' and r['eligible']=='1' and r['named']=='1'}
  for f in sorted((WORK/'12S/downloads').glob('*.gz')):
   with gzip.open(f,'rt') as h:
    for r in SeqIO.parse(h,'genbank'):
     if r.id not in audit:continue
     row=audit[r.id];feat=next((x for x in r.features if x.type==row['method'] and str(x.location)==row['location']),None);s=str(feat.extract(r.seq) if feat else r.seq).upper();assert len(s)==int(row['length']);seq[r.id]=s;meta.append(row)
 with (out/'panel.csv').open('w') as h:
  w=csv.DictWriter(h,meta[0]);w.writeheader();w.writerows(meta)
 (out/'target.fasta').write_text(''.join('>'+k+'\n'+v+'\n' for k,v in seq.items()))
 with (out/'aligned.fasta').open('w') as h:subprocess.run(['mafft','--quiet','--adjustdirection','--retree','2','--maxiterate','0','--thread','4',str(out/'target.fasta')],stdout=h,check=True)
 ps=rows(ROOT/gene/'candidate_results.csv');(out/'independent_primers.tsv').write_text(''.join('\t'.join([p['name'],p['F'],p['R']])+'\n' for p in ps));(out/'background.fasta').write_text('')
 print(gene,len(seq),'records',len({r['organism'] for r in meta}),'species',flush=True)
