from pathlib import Path
from Bio import SeqIO
import csv,json,itertools,collections
r=Path(__file__).resolve().parent
read=lambda p:list(csv.DictReader(p.open()))
def save(fn,rows):
 if not rows:return
 with (r/fn).open('w') as h:w=csv.DictWriter(h,rows[0]);w.writeheader();w.writerows(rows)
ps=read(r/'candidates.csv');meta={x['accession']:x for x in read(r/'panel.csv')};al={x.id.removeprefix('_R_'):str(x.seq).upper() for x in SeqIO.parse(r/'target.aligned.fasta','fasta')};ref=json.loads((r/'reference.json').read_text())['accession'];co=[i for i,b in enumerate(al[ref]) if b!='-']
hits=list(csv.DictReader((r/'independent_hits.tsv').open(),delimiter='\t'));hit={(x['name'],x['role'],x['model'],x['accession']):x for x in hits};out=[];details=[];bars={};families=[]
for p in ps:
 name=p['name'];bar={};eligible=[]
 for acc,s in al.items():
  if name=='Anth28S' or not p['F_start']:continue
  a,b,c,d=[co[int(p[k])-1] for k in ['F_start','F_end','R_start','R_end']];f=s[a:b+1].replace('-','');z=s[c:d+1].replace('-','');inner=s[b+1:c].replace('-','');good=len(f)==len(p['F']) and len(z)==len(p['R']) and set(f+z)<=set('ACGT') and 50<=len(inner) and 100<=len(f+inner+z)<=500
  if good:
   eligible.append(acc)
   if set(inner)<=set('ACGT'):bar[meta[acc]['species']]=inner
  h=hit.get((name,'target','WC',acc),{});details.append(dict(name=name,accession=acc,species=meta[acc]['species'],family=meta[acc]['family'],evaluable=int(good),covered=int(good and h.get('le2last3')=='1'),exact=int(good and h.get('exact')=='1')))
 if name=='Anth28S' or not p['F_start']:
  eligible=list(al)
  bb=collections.defaultdict(set)
  for v in csv.DictReader((r/'target_barcodes.tsv').open(),delimiter='\t'):
   if v['name']==name:bb[v['accession']].add(v['sequence'])
  bar={meta[a]['species']:next(iter(v)) for a,v in bb.items() if len(v)==1}

 bars[name]=bar;cov=sum(hit.get((name,'target','WC',a),{}).get('le2last3')=='1' for a in eligible);exact=sum(hit.get((name,'target','WC',a),{}).get('exact')=='1' for a in eligible)
 x=dict(p,n=len(eligible),covered=cov,exact=exact,barcode_n=len(bar),barcode_unique=sum(v not in [v2 for k2,v2 in bar.items() if k2!=k] for k,v in bar.items()))
 for model,key,label in [('WC','le2last3','WC'),('GT','le2last3','GT'),('WC','le3','WC3'),('GT','le2','GT2')]:x['bg_'+label]=sum(h[key]=='1' for h in hits if h['name']==name and h['role']=='background' and h['model']==model)
 out.append(x)
 for group in sorted({m['group'] for m in meta.values() if m['group']!='Actiniaria'}):
  mm=[a for a,m in meta.items() if m['group']==group];families.append(dict(name=name,group=group,n=len(mm),WC=sum(hit.get((name,'background','WC',a),{}).get('le2last3')=='1' for a in mm),GT=sum(hit.get((name,'background','GT',a),{}).get('le2last3')=='1' for a in mm)))
save('comparison.csv',out);save('species_results.csv',details);save('background_by_group.csv',families)
familyrows=[]
for name in bars:
 for family in sorted({x['family'] for x in details}):
  dd=[x for x in details if x['name']==name and x['family']==family]
  if dd:familyrows.append(dict(name=name,family=family,total=len(dd),evaluable=sum(x['evaluable'] for x in dd),covered=sum(x['covered'] for x in dd),exact=sum(x['exact'] for x in dd)))
save('target_by_family.csv',familyrows)
# Every mitochondrial same-barcode species pair is tested against each 28S window.
covered={(x['name'],x['species']):x['covered'] for x in details}
old=r/'mitochondrial_baseline';pairs=[];mitobars={}
for gene in ['12S','16S']:
 d=old/gene;aa={x.id.removeprefix('_R_'):str(x.seq).upper() for x in SeqIO.parse(d/'aligned.fasta','fasta')};cc=[i for i,b in enumerate(aa['NC_000933.1']) if b!='-'];mm={x['accession']:x['organism'] for x in read(d/'target_panel.csv')}
 elig=collections.defaultdict(set)
 for e in read(d/'species_results.csv'):
  if e['evaluable']=='1':elig[e['name']].add(e['accession'])
 for p in read(d/'candidate_results.csv'):
  b,c=cc[int(p['F_end'])-1],cc[int(p['R_start'])-1];groups=collections.defaultdict(set);mb={}
  for acc,s in aa.items():
   t=s[b+1:c].replace('-','')
   if acc in mm and acc in elig[p['name']] and t and set(t)<=set('ACGT'):groups[t].add(mm[acc]);mb[mm[acc]]=t
  mitobars[gene,p['name']]=mb
  for gg in groups.values():
   for a,b in itertools.combinations(sorted(gg),2):
    for name,bar in bars.items():
     status='missing_28S' if a not in bar or b not in bar else 'resolved' if bar[a]!=bar[b] else 'still_shared';pairs.append(dict(gene=gene,mitochondrial_pair=p['name'],species_a=a,species_b=b,candidate_28S=name,status=status,both_covered=int(bool(covered.get((name,a))) and bool(covered.get((name,b)))) ))
classification={}
for v in pairs:
 key=v['species_a'],v['species_b']
 if key not in classification:
  a,b=key;compared=[(g,bb[a]==bb[b]) for (g,n),bb in mitobars.items() if a in bb and b in bb]
  classification[key]=int({g for g,eq in compared}=={'12S','16S'} and all(eq for g,eq in compared))
 v['all_available_12S_16S_windows_shared']=classification[key]
save('mitochondrial_conflict_rescue.csv',pairs)
core={(v['candidate_28S'],v['species_a'],v['species_b']):v for v in pairs if v['all_available_12S_16S_windows_shared']}
save('core_conflict_rescue.csv',list(core.values()))
summ=[]
for x in out:
 rr=[v for v in pairs if v['candidate_28S']==x['name']];uniq={(v['species_a'],v['species_b']):v['status'] for v in rr};summ.append(dict(name=x['name'],mitochondrial_conflict_pairs=len(uniq),evaluable=sum(v!='missing_28S' for v in uniq.values()),resolved=sum(v=='resolved' for v in uniq.values()),still_shared=sum(v=='still_shared' for v in uniq.values()),unresolved_across_both_markers=sum(classification[k] for k in uniq),core_resolved_and_covered=len({(v['species_a'],v['species_b']) for v in rr if v['all_available_12S_16S_windows_shared'] and v['status']=='resolved' and v['both_covered']}),resolved_and_both_covered=len({(v['species_a'],v['species_b']) for v in rr if v['status']=='resolved' and v['both_covered']})))
save('rescue_summary.csv',summ)
(r/'barcodes.json').write_text(json.dumps(bars,indent=2)+'\n')
print('summarized',len(out),'candidates',flush=True)
