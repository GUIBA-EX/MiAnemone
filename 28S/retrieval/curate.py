from pathlib import Path
import xml.etree.ElementTree as E,re,csv,json
from collections import defaultdict
from Bio.Seq import Seq
r=Path(__file__).resolve().parent;rows=[];seqs={};excluded=[]
for f in sorted((r/'data').glob('*_0*.xml')):
 group=f.name.split('_')[0]
 for x in E.parse(f).getroot().findall('GBSeq'):
  acc=x.findtext('GBSeq_accession-version');sp=x.findtext('GBSeq_organism');raw=x.findtext('GBSeq_sequence','').upper();tax=x.findtext('GBSeq_taxonomy','');definition=x.findtext('GBSeq_definition','');features=[]
  for ft in x.findall('./GBSeq_feature-table/GBFeature'):
   qq={q.findtext('GBQualifier_name'):q.findtext('GBQualifier_value','') for q in ft.findall('./GBFeature_quals/GBQualifier')}
   if ft.findtext('GBFeature_key')=='rRNA' and re.search(r'28S|large subunit',qq.get('product',''),re.I):
    parts=[]
    for iv in ft.findall('./GBFeature_intervals/GBInterval'):
     a=int(iv.findtext('GBInterval_from',iv.findtext('GBInterval_point','0')));b=int(iv.findtext('GBInterval_to',str(a)))
     if a and b:parts.append(raw[min(a,b)-1:max(a,b)])
    s=''.join(parts)
    if 'complement' in ft.findtext('GBFeature_location',''):s=str(Seq(s).reverse_complement())
    if s:features.append(s)
  if features:s=max(features,key=len);source='annotated_rRNA'
  elif x.findtext('GBSeq_moltype') not in ['mRNA','protein'] and re.search(r'28S.*ribosomal (?:RNA|DNA|gene)|large subunit ribosomal (?:RNA|DNA|gene)',definition,re.I) and not re.search(r'18S|5\.8S|internal transcribed',definition,re.I):s=raw;source='single_locus_title'
  else:excluded.append([acc,group,sp,'no_unambiguous_28S_feature']);continue
  if not re.fullmatch('[A-Z][a-z-]+ [a-z][a-z-]+',sp):excluded.append([acc,group,sp,'non_binomial']);continue
  if len(s)<150 or sum(c not in 'ACGT' for c in s)/len(s)>.02:excluded.append([acc,group,sp,'short_or_ambiguous']);continue
  family=next((v.strip() for v in tax.split(';') if v.strip().endswith('idae')),'unassigned')
  rows.append(dict(accession=acc,species=sp,group=group,family=family,length=len(s),concrete=sum(c in 'ACGT' for c in s),extraction=source,taxonomy=tax));seqs[acc]=s
by=defaultdict(list)
for x in rows:by[x['group'],x['species']].append(x)
selected=[max(v,key=lambda x:(x['concrete'],x['length'],x['accession'])) for v in by.values()]
def csvout(fn,rr):
 with (r/fn).open('w') as h:w=csv.DictWriter(h,rr[0]);w.writeheader();w.writerows(rr)
csvout('all_records.csv',rows);csvout('panel.csv',selected)
with (r/'exclusions.csv').open('w') as h:w=csv.writer(h);w.writerow(['accession','group','species','reason']);w.writerows(excluded)
for role in ['target','background']:
 rr=[x for x in selected if (x['group']=='Actiniaria')==(role=='target')];(r/(role+'.fasta')).write_text(''.join('>'+x['accession']+' '+x['species']+'\n'+seqs[x['accession']]+'\n' for x in rr))
target=[x for x in selected if x['group']=='Actiniaria'];ref=max(target,key=lambda x:x['concrete']);(r/'reference.fasta').write_text('>'+ref['accession']+'\n'+seqs[ref['accession']]+'\n');(r/'reference.json').write_text(json.dumps(ref,indent=2)+'\n');print('records',len(rows),'representatives',len(selected),'target',len(target),'reference',ref,flush=True)
