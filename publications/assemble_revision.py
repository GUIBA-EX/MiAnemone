"""Append evidence tables to the revision body; no manual transcription of metrics."""
import csv
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def rows(p):return list(csv.DictReader((ROOT/p).open()))
def interval(a,b,d):
 a,b=float(a),float(b)
 return f'{a:.{d}f}' if a==b else f'{a:.{d}f}–{b:.{d}f}'
def table(headers,rr):return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(map(str,r))+' |' for r in rr])+'\n'
s=(ROOT/'publications/revision_body.md').read_text();properties=rows('revision/oligo_properties.csv');rr=[]
for r in properties:
 if r['scenario']!='50_per_variant':continue
 seq=r['sequence']
 if (r['gene'],r['pair'],r['side'])==('12S','XS','F'):seq=seq[:9]+'**'+seq[9]+'**'+seq[10:]
 rr.append([r['gene'],r['pair']+'-'+r['side'],seq,r['length'],r['variants'],interval(r['GC_min'],r['GC_max'],0),interval(r['Tm_min'],r['Tm_max'],1)])
s+='\n## Table 1. Oligos and predicted properties\n\nAll sequences run 5′–3′. The adopted G10A base is bold. Tm is at 50 nM per concrete variant; the 200 nM-per-direction scenario is supplied in `revision/oligo_properties.csv`. Variant count is per end, not the product of forward and reverse degeneracies.\n\n'+table(['Gene','Oligo','Sequence','nt','Variants','GC %','Tm °C'],rr)
common={(r['gene'],r['pair']):r for r in rows('revision/panel_comparison.csv') if r['panel']=='common'};length={(r['gene'],r['pair']):r for r in rows('revision/length_summary.csv')};joint={(r['gene'],r['pair']):r for r in rows('revision/joint_resolution.csv') if r['panel']=='common'}
rr=[]
for p in rows('all_results.csv'):
 k=p['gene'],p['name'];c=common[k];l=length[k];j=joint[k];rr.append([*k,f"{p['covered']}/{p['evaluable']} ({100*int(p['covered'])/int(p['evaluable']):.1f}%)",f"{p['exact']}/{p['evaluable']}",f"{c['covered']}/{c['n']}",f"{c['exact']}/{c['n']}",f"{l['min']}/{float(l['median']):g}/{l['max']}",f"{j['covered_unique']}/{j['n']}"])
s+='\n## Table 2. Coverage, exact matching, product lengths and joint discrimination\n\nOwn denotes the window-specific species panel; common denotes the same evaluable species across all windows. Product lengths are minimum/median/maximum among distinct accession–length predictions under strict WC. Joint uniqueness uses the smaller common panel with complete inserts at every window. Computational priorities are discussed above; no pair is wet-laboratory validated.\n\n'+table(['Gene','Pair','Own covered (%)','Own exact','Common covered','Common exact','Product bp','Covered + unique'],rr)
s+='\n## Table 3. Factorial decomposition of non-nested pairs\n\nFixed changes and degeneracy changes are applied separately (1 = applied, 0 = absent). The 0/1 state is the nested control; the 1/1 state is the preceding D candidate. The 12S XS 0/1 state is now the released XSD; the other two D pairs retain their 1/1 states. Thermal pass requires minimum variant Tm ≥45°C, structure Tm <55°C, and maximum F/R difference ≤5°C in both concentration scenarios. The adopted 12S XSD is prioritized for short-fragment evaluation.\n\n'+table(['Gene','D pair','Fixed','Degenerate','Covered','Exact','Evaluable','Terminal','Thermal'],[[r[k] for k in ['gene','pair','fixed_changes','degenerate_changes','covered','exact','n']]+['Pass' if r[k]=='True' else 'Fail' for k in ['terminal_pass','thermal_pass']] for r in rows('revision/prior_XSD/factorial_filters.csv')])
support={(r['gene'],r['pair'],r['edits']):r for r in rows('revision/threshold_sensitivity.csv')};rr=[]
for r in rows('audit/consensus_trials.csv'):
 t=support[r['gene'],r['pair'],r['edits']];rr.append([r['gene'],r['pair'],r['edits'],f"{100*float(t['minimum_frequency']):.1f}",r['covered']+'/'+r['exact'],r['lost_covered'],r['added_background_model_hits'],f"{float(r['max_delta_Tm']):.2f}"])
s+='\n## Table 4. All 23 post-G10A consensus alternatives\n\nThis audit predates the 12S XSD replacement; XSD edits refer to the preceding pair. Frequency is the minimum new-base frequency across the edited positions. Lost counts track identities rather than net change. Added background counts sum added accession–model hits across the four original models and can count one accession more than once. The supplement gives all threshold outcomes, terminal and structure flags, and accession-level gains/losses.\n\n'+table(['Gene','Pair','Edits','Freq %','Covered/exact','Lost','Added bg','Max ΔTm'],rr)
main,audit=s.split('## Table 4.',1)
(ROOT/'publications/MiAnemone_supplement.md').write_text('# MiAnemone supplementary material\n\n**Fei Xia**\n\n## Table S1.'+audit)
(ROOT/'publications/MiAnemone_preprint.md').write_text(main+'\n[Supplementary candidate audit](MiAnemone_supplement.md)\n')
