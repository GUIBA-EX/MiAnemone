import csv,json,math
from pathlib import Path
from collections import defaultdict
from verify import rows,fasta
r=Path(__file__).resolve().parents[1];checks=[]
for gene in ['12S','16S']:
 hits={(x['name'],x['role'],x['model'],x['accession']):int(x['le2last3']) for x in csv.DictReader((r/'revision'/('independent_'+gene)/'independent_hits.tsv').open(),delimiter='\t')}
 rr=[x for x in rows(r/'revision/unified_record_hits.csv') if x['gene']==gene and x['max_mismatches']=='2'];bad=[x for x in rr if int(x['hit'])!=hits.get((x['pair'],x['role'],x['model'],x['accession']),0)]
 assert not bad,bad[:3]
 checks.append(dict(gene=gene,record_model_checks=len(rr),disagreements=len(bad)))
(r/'revision/independent_validation.json').write_text(json.dumps(checks,indent=2)+'\n')
scores=defaultdict(float)
for x in rows(r/'revision/mismatch_positions.csv'):
 scores[x['gene'],x['pair'],x['accession']]+=int(x['mismatch'])*(1+4*math.exp(-int(x['distance_from_3prime'])/3))
with (r/'revision/positional_penalties.csv').open('w') as f:
 w=csv.writer(f);w.writerow(['gene','pair','accession','positional_penalty']);w.writerows([*k,v] for k,v in scores.items())
print(checks)
