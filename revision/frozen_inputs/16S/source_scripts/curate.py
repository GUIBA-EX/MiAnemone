"""Extract strand-aware mitochondrial rrnL features and create species representatives."""
from pathlib import Path
import csv,json,gzip,re,collections,hashlib
from Bio import SeqIO
p=Path(__file__).resolve().parents[1];audit=[];seqs={};source_meta={}
def tab(path,rr):
 with path.open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rr[0]));w.writeheader();w.writerows(rr)
def fasta(path,ids):
 with path.open('w') as f:
  for a in ids:f.write('>'+a+'\n'+seqs[a]+'\n')
def locus(s,mito):
 s=s.lower()
 if re.search(r'\b(18s|28s|12s|protein|methyltransferase)\b',s):return False
 return bool(re.search(r'\b(16s|rrnl|rnl|l-rrna)\b',s) or (mito and 'large subunit' in s))
for f in sorted((p/'downloads').glob('*.gb.gz')):
 with gzip.open(f,'rt') as stream:
  for r in SeqIO.parse(stream,'genbank'):
   tax=r.annotations.get('taxonomy',[]);org=r.annotations.get('organism','');desc=r.description;source=next((x.qualifiers for x in r.features if x.type=='source'),{});mito='mitochond' in (desc+str(source.get('organelle',''))).lower();tx=next((s[6:] for s in source.get('db_xref',[]) if s.startswith('taxon:')),'');matches=[];status='';method='';loc='';seq=''
   if 'Actiniaria' not in tax:status='outside_target_lineage'
   elif not r.seq.defined:status='undefined_sequence'
   else:
    for x in r.features:
     if x.type not in ['rRNA','gene'] or not x.location:continue
     label=' '.join(x.qualifiers.get('gene',[])+x.qualifiers.get('product',[])+x.qualifiers.get('standard_name',[]))
     if locus(label,mito):matches.append((x.type=='rRNA',len(x),x))
    if matches:
     x=max(matches,key=lambda t:t[:2])[2];seq=str(x.extract(r.seq)).upper();loc=str(x.location);method='feature_'+x.type
    elif locus(desc,mito) and len(r)<=3000 and not re.search('tRNA|transfer RNA|cytochrome|subunit I |ND[1-6]|complete genome|internal transcribed| and ',desc,re.I) and all(x.type not in ('CDS','tRNA') for x in r.features):
     seq=str(r.seq).upper();loc='whole_record';method='single_locus_title'
    else:status='no_supported_16S_feature'
   amb=sum(c not in 'ACGT' for c in seq)
   if seq:
    status='eligible' if len(seq)>=120 and amb<=0.05*len(seq) and set(seq)<=set('ACGTRYSWKMBDHVN') else 'length_or_ambiguity'
   named=bool(re.fullmatch(r'[A-Z][a-z-]+ (?:\([A-Z][a-z-]+\) )?[a-z][a-z-]+',org))
   row=dict(accession=r.id,organism=org,taxid=tx,family=next((v for v in reversed(tax) if v.endswith('idae')),''),lineage=';'.join(tax),named_binomial=int(named),method=method,mitochondrial_source_explicit=int(mito),location=loc,length=len(seq),ambiguous=amb,status=status,representative=0,description=desc)
   audit.append(row)
   if status=='eligible':seqs[r.id]=seq;source_meta[r.id]=row
groups=collections.defaultdict(list)
for a,x in source_meta.items():
 if x['named_binomial']:groups[x['organism']].append(a)
reps=[]
for name,ids in sorted(groups.items()):
 a=max(ids,key=lambda a:(source_meta[a]['length']-source_meta[a]['ambiguous'],source_meta[a]['method'].startswith('feature'),-source_meta[a]['ambiguous'],a.startswith('NC_'),a));source_meta[a]['representative']=1;reps.append(a)
assert len({x['accession'] for x in audit})==len(audit)==json.loads((p/'downloads/search.json').read_text())['count']
fasta(p/'all_eligible.fasta',sorted(seqs));fasta(p/'representatives.fasta',reps);tab(p/'accession_audit.csv',audit);tab(p/'species_panel.csv',[source_meta[a] for a in reps])
summary=dict(downloaded=len(audit),eligible=len(seqs),named_species=len(reps),genera=len({s.split()[0] for s in groups}),families=len({source_meta[a]['family'] for a in reps if source_meta[a]['family']}),status_counts=collections.Counter(x['status'] for x in audit),method_counts=collections.Counter(x['method'] for x in source_meta.values()))
(p/'summary.json').write_text(json.dumps(summary,indent=2));print(summary)
