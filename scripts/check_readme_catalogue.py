"""Validate the integrated three-marker README catalogue and binding maps."""
from pathlib import Path
import csv,xml.etree.ElementTree as ET
r=Path(__file__).resolve().parents[1]
rows=list(csv.DictReader((r/'all_results.csv').open()))
for gene in ('12S','16S','28S'):
 tree=ET.parse(r/f'figures/{gene}_binding.svg');text=' '.join(tree.getroot().itertext())
 assert ('MT676790.1' if gene=='28S' else 'NC_000933.1') in text
 for p in [x for x in rows if x['gene']==gene]:
  assert f"F {p['F_start']}–{p['F_end']}" in text
  assert f"R {p['R_start']}–{p['R_end']}" in text
  assert f"{p['span']} bp" in text
for file in ('README.md','README.en.md'):
 s=(r/file).read_text()
 for gene in ('12S','16S','28S'):assert f'figures/{gene}_binding.svg' in s
 for p in rows:
  lines=[line for line in s.splitlines() if line.startswith('| '+p['name']+' |') and p['F'] in line and p['R'] in line]
  assert len(lines)==1,(file,p['gene'],p['name'])
  assert f"figures/template_matching/{p['gene']}_{p['name']}.svg" in lines[0]
  for k in ('covered','exact'):assert f"{p[k]}/{p['evaluable']}" in lines[0]
 assert s.index('### 12S\n')<s.index('### 16S\n')<s.index('### 28S\n')
 assert '## Nuclear 28S / 核28S' not in s
 assert s.count('<details>')==s.count('</details>')
print('PASS: bilingual 30-pair catalogue, per-pair diagram links and three reference-coordinate maps')
