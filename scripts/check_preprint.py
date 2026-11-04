"""Validate the revision's sequence, evidence tables and local Markdown targets."""
import csv,re,json,sys,hashlib
from pathlib import Path
root=Path(__file__).resolve().parents[1]
def rows(p):return list(csv.DictReader((root/p).open()))
pairs=rows('all_results.csv');paper=(root/'publications/MiAnemone_preprint.md').read_text();plain=paper.replace('**','');fa=(root/'primers.fasta').read_text().splitlines();expected={f"{r['gene']}_{r['name']}_{e}":r[e] for r in pairs for e in 'FR'}
assert dict(zip([k[1:] for k in fa[::2]],fa[1::2],strict=True))==expected
assert len(pairs)==22 and len(expected)==44
assert all(q in plain for q in expected.values())
assert 'universal PCR primers' not in paper.splitlines()[0]
assert 'Supplementary expanded-background analyses' not in paper
assert len(rows('revision/prior_XSD/factorial_controls.csv'))==12
assert len(rows('audit/consensus_trials.csv'))==23
assert len(rows('revision/threshold_sensitivity.csv'))==368
checks=json.loads((root/'revision/independent_validation.json').read_text());assert sum(x['record_model_checks'] for x in checks)==57584 and not any(x['disagreements'] for x in checks)
for p in [root/'README.md',root/'README.en.md',root/'publications/MiAnemone_preprint.md',root/'publications/response_to_review.md']:
 for url in re.findall(r'\]\(([^)]+)\)',p.read_text()):
  if '://' in url or url.startswith('#'):continue
  assert (p.parent/url.split('#')[0]).exists(),(p,url)
for g in ['12S','16S']:
 for r in rows(g+'/candidate_results.csv'):
  assert (root/f"figures/template_matching/{g}_{r['name']}.svg").exists()
if len(sys.argv)>1:
 text=Path(sys.argv[1]).read_text();assert all(s in text for s in expected.values()),'Missing sequence in rendered PDF'
summary=dict(pair_count=22,oligo_count=44,factorial_controls=12,substitution_candidates=23,threshold_combinations=368,strict_independent_checks=57584,disagreements=0,midori_included=False,status='revision draft; outstanding items documented in response')
(root/'revision/document_checks.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary)
