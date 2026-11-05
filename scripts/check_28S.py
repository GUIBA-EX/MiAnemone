from pathlib import Path
import csv,json,math,xml.etree.ElementTree as ET
from Bio.Data.IUPACData import ambiguous_dna_values as I
r=Path(__file__).resolve().parents[1];d=r/'28S';read=lambda f:list(csv.DictReader(f.open()));pp=read(d/'candidate_results.csv');ss=read(d/'selected_species_results.csv');assert {x['name'] for x in pp}=={'XS','XSD','S','SD','M','MD','L','LD'}
for p in pp:
 rr=[x for x in ss if x['name']==p['name']]
 for k in ['evaluable','covered','exact']:assert sum(int(x[k]) for x in rr)==int(p[k]),(p['name'],k)
 for side in 'FR':assert len(p[side])==int(p[side+'_end'])-int(p[side+'_start'])+1 and math.prod(len(I[b]) for b in p[side])==int(p['deg_'+side])<=8
 if p['name'].endswith('D'):
  a=next(x for x in pp if x['name']==p['name'][:-1]);assert all(a[k]==p[k] for k in ['F_start','F_end','R_start','R_end'])
 for suffix in ['', '_preview']:
  f=r/'figures/template_matching'/f"28S_{p['name']}{suffix}.svg";ET.parse(f)
  assert p['F'] in f.read_text() and p['R'] in f.read_text()
 for file in [r/'README.md',r/'README.en.md',r/'publications/MiAnemone_preprint.md']:
  s=file.read_text();assert p['F'] in s and p['R'] in s;assert f"{p['covered']}/{p['evaluable']}" in s and f"{p['exact']}/{p['evaluable']}" in s
v=json.loads((d/'verification.json').read_text());assert v['checks']==90984 and not v['disagreements']
print('PASS 28S: eight pairs, fixed coordinates, degeneracy, species totals, 16 diagrams, manuscript sequences and independent checks')
