"""Recompute plotted counts and entropy from accession-level inputs."""
import csv,math,json
from pathlib import Path
from collections import Counter
from verify import rows,fasta
r=Path(__file__).resolve().parents[1];report={}
for g in ['12S','16S']:
 ps=rows(r/g/'candidate_results.csv');detail=rows(r/g/'species_results.csv');al={k.removeprefix('_R_'):s for k,s in fasta(r/g/'aligned.fasta').items()};meta=rows(r/g/'target_panel.csv');co=[i for i,c in enumerate(al['NC_000933.1']) if c!='-'];er=[x for x in rows(r/'revision/entropy.csv') if x['gene']==g]
 for j,x in enumerate(er):
  ct=Counter(al[z['accession']][co[j]] for z in meta if al[z['accession']][co[j]] in 'ACGT');n=sum(ct.values());h=-sum(v/n*math.log2(v/n) for v in ct.values()) if n else 0
  assert n==int(x['n']) and abs(h-float(x['entropy']))<1e-10
 by={p['name']:{x['accession']:x for x in detail if x['name']==p['name']} for p in ps};ids=set.intersection(*[{k for k,x in d.items() if x['evaluable']=='1'} for d in by.values()]);complete=set(ids)
 for p in ps:
  b,c=co[int(p['F_end'])-1],co[int(p['R_start'])-1];complete&={k for k in ids if (q:=al[k][b+1:c].replace('-','')) and set(q)<=set('ACGT')}
 for x in rows(r/'revision/panel_comparison.csv'):
  if x['gene']!=g or x['panel']!='common':continue
  assert int(x['n'])==len(ids)
  for key in ['covered','exact']:assert int(x[key])==sum(int(by[x['pair']][k][key]) for k in ids)
 for x in rows(r/'revision/joint_resolution.csv'):
  if x['gene']!=g or x['panel']!='common':continue
  p=next(p for p in ps if p['name']==x['pair']);b,c=co[int(p['F_end'])-1],co[int(p['R_start'])-1];seq={k:al[k][b+1:c].replace('-','') for k in complete};cnt=Counter(seq.values());assert int(x['covered_unique'])==sum(cnt[v]==1 and by[x['pair']][k]['covered']=='1' for k,v in seq.items())
 groups=['Scleractinia','Antipatharia','Corallimorpharia','Zoantharia','Octocorallia','Hydrozoa'];hits=rows(r/g/'background_hits.csv')
 for x in hits:assert sum(z in x['lineage'].split(';') for z in groups)==1
 for p in ps:
  for m,a,b in [('WC','WC','WC3'),('GT','GT','GT2')]:
   hh=[x for x in hits if x['name']==p['name'] and x['model']==m];assert sum(int(x['strict']) for x in hh)==int(p[a]);assert sum(int(x['relaxed']) for x in hh)==int(p[b])
 for p in ps:
  if not p['name'].endswith('D'):
   d=next(q for q in ps if q['name']==p['name']+'D');assert all(p[k]==d[k] for k in ['F_start','F_end','R_start','R_end'])
 report[g]=dict(entropy_positions=len(er),called_species_min=min(int(x['n']) for x in er),called_species_max=max(int(x['n']) for x in er),common_coverage_species=len(ids),common_insert_species=len(complete),pairs=len(ps),background_hits_accounted=True)
from Bio.Data.IUPACData import ambiguous_dna_values as I
pairs={(p['gene'],p['name']):p for p in rows(r/'all_results.csv')}
mm=rows(r/'revision/mismatch_positions.csv')
for x in mm:
 q=pairs[x['gene'],x['pair']][x['side']];j=len(q)-1-int(x['distance_from_3prime']);assert q[j]==x['primer'];assert int(x['mismatch'])==int(x['template'] not in I[q[j]])
curves={}
for x in rows(r/'revision/unified_models.csv'):
 if x['terminal_zero']=='3':curves.setdefault((x['gene'],x['pair'],x['role'],x['model']),{})[int(x['max_mismatches'])]=int(x['matched'])
for c in curves.values():assert set(c)==set(range(5)) and all(c[i]<=c[i+1] for i in range(4))
report['mismatch_records_checked']=len(mm);report['monotonic_model_curves']=len(curves)
for name in ['review_common_panel','review_entropy','review_mismatch_profiles','review_unified_models','review_background_all','review_joint_resolution']:
 s=(r/'publications/figures'/f'{name}.svg').read_text();assert 'Avenir Next' in s
for path in [r/'README.md',r/'README.en.md',r/'publications/MiAnemone_preprint.md']:
 s=path.read_text();assert 'MIDORI' not in s
 for name in ['review_common_panel','review_entropy','review_mismatch_profiles','review_unified_models','review_background_all','review_joint_resolution']:assert name+'.svg' in s
(r/'revision/figure_checks.json').write_text(json.dumps(report,indent=2)+'\n');print(report)
