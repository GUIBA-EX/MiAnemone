from pathlib import Path
import csv,gzip,collections,re,json
from Bio import SeqIO
p=Path(__file__).resolve().parents[1];d=p/'12S';old=p.parents[1]/'work/scleractinia16s/downloads';seen=set();seq={};records=[];groups=collections.defaultdict(list)
def write(path,rows):
 with path.open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
for f,role in [(x,'target') for x in sorted((p/'downloads').glob('*.gz'))]+[(x,'background') for x in sorted(old.glob('*.gb'))]:
 with (gzip.open(f,'rt') if f.suffix=='.gz' else f.open()) as h:
  for r in SeqIO.parse(h,'genbank'):
   if (role,r.id) in seen:continue
   seen.add((role,r.id));tax=r.annotations.get('taxonomy',[]);org=r.annotations.get('organism','');desc=r.description.lower();src=next((x.qualifiers for x in r.features if x.type=='source'),{});mito='mitochond' in desc+str(src.get('organelle','')).lower()
   if role=='background' and ('Cnidaria' not in tax or 'Actiniaria' in tax):continue
   if role=='target' and 'Actiniaria' not in tax:continue
   s='';method='';loc='';fs=[]
   if r.seq.defined and not r.id.startswith(('XM_','XR_','NM_')) and not re.search('methyltransferase|pseudogene',desc):
    for feat in r.features:
     if feat.type not in ('rRNA','gene') or not feat.location:continue
     t=' '.join(feat.qualifiers.get('gene',[])+feat.qualifiers.get('product',[])).lower()
     if re.search(r'\b(18s|16s|protein)\b',t):continue
     if '12s' in t or (mito and re.search(r'\b(rrns|rns|s-rrna)\b|small subunit',t)):fs.append(feat)
    if fs:
     feat=max(fs,key=lambda x:(x.type=='rRNA',len(x)));s=str(feat.extract(r.seq)).upper();method=feat.type;loc=str(feat.location)
    elif '12s' in desc and len(r)<=1600 and not re.search(r'trna|16s|18s|complete genome|protein|mrna| and ',desc) and not any(x.type in ('tRNA','CDS') for x in r.features):s=str(r.seq).upper();method='title';loc='whole'
   ok=120<=len(s)<=2000 and sum(b not in 'ACGT' for b in s)<=.05*len(s);named=bool(re.fullmatch('[A-Z][a-z-]+ (?:\\([A-Z][a-z-]+\\) )?[a-z][a-z-]+',org))
   row=dict(accession=r.id,organism=org,family=next((x for x in reversed(tax) if x.endswith('idae')),''),lineage=';'.join(tax),role=role,length=len(s),location=loc,method=method,eligible=int(ok),named=int(named))
   records.append(row)
   if ok:
    seq[r.id]=s
    if named:groups[role,org].append(row)
for role in ['target','background']:
 rr=[max(v,key=lambda r:('complete genome' in r.get('description','') and r['method']=='rRNA',sum(c in 'ACGT' for c in seq[r['accession']]),r['method']=='rRNA',r['accession'].startswith('NC_'),r['accession'])) for (ro,org),v in groups.items() if ro==role]
 write(d/(role+'_panel.csv'),rr);(d/(role+'.fasta')).write_text(''.join('>'+r['accession']+'\n'+seq[r['accession']]+'\n' for r in rr));print(role,len(rr))
 if role=='target':
  anchor='NC_000933.1';assert anchor in seq;(d/'anchor.txt').write_text(anchor+'\n');(d/'align_input.fasta').write_text('>'+anchor+'\n'+seq[anchor]+'\n'+''.join('>'+r['accession']+'\n'+seq[r['accession']]+'\n' for r in rr if r['accession']!=anchor))
write(p/'audit.csv',records)
