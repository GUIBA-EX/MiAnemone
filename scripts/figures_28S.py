from pathlib import Path
import json
from template_matching import render,ROOT,rows,fasta
p=ROOT/'28S';al={k.removeprefix('_R_'):v for k,v in fasta(p/'target.aligned.fasta').items()};ref=json.loads((p/'reference.json').read_text())['accession'];cols=[i for i,b in enumerate(al[ref]) if b!='-'];rr=rows(p/'selected_species_results.csv')
for pair in rows(p/'candidate_results.csv'):
 eligible=[x for x in rr if x['name']==pair['name'] and x['evaluable']=='1'];assert len(eligible)==int(pair['evaluable'])
 render('28S',pair,eligible,al,cols)
 render('28S',pair,eligible,al,cols,preview=True)
