"""Reproducible review diagnostics. No replacement of released oligos."""
import csv,math,re,itertools,statistics
from collections import Counter,defaultdict
from pathlib import Path
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
from Bio.Data.IUPACData import ambiguous_dna_values as I
from verify import rows,fasta
from independent import flags
import primer3
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'revision'
def write(n,rr):
 with (OUT/(n+'.csv')).open('w') as f:
  w=csv.DictWriter(f,rr[0],lineterminator='\n');w.writeheader();w.writerows(rr)
def expand(s):return list(map(''.join,itertools.product(*(I[c] for c in s))))
kw=dict(mv_conc=50,dv_conc=1.5,dntp_conc=.6)
res=defaultdict(list);oligos=[]
aligner=PairwiseAligner(mode='global',match_score=0,mismatch_score=-1,open_gap_score=-1,extend_gap_score=-1)
for gene in ['12S','16S']:
 d=ROOT/gene;ps=rows(d/'candidate_results.csv');meta={r['accession']:r for r in rows(d/'target_panel.csv')};al={k.removeprefix('_R_'):s for k,s in fasta(d/'aligned.fasta').items()};cols=[i for i,c in enumerate(al['NC_000933.1']) if c!='-'];detail=rows(d/'species_results.csv');raw=fasta(d/'target.fasta')
 for p in ps:
  name=p['name'];dd={r['accession']:r for r in detail if r['name']==name};ins={}
  a,b,c,e=[cols[int(p[k])-1] for k in ['F_start','F_end','R_start','R_end']]
  for acc in meta:
   s=al[acc];spans=bool(s[:a+1].strip('-') and s[e:].strip('-'));sites=s[a:b+1]+s[c:e+1];project=''.join(s[j] for j in cols[int(p['F_start'])-1:int(p['F_end'])]+cols[int(p['R_start'])-1:int(p['R_end'])]);reason='not_spanning' if not spans else 'ambiguous_binding' if set(sites)-set('ACGT-') else 'evaluable'
   res['template_status'].append(dict(gene=gene,pair=name,accession=acc,species=meta[acc]['organism'],family=meta[acc]['family'],status=reason,projected_gap=('-' in project),aligned_site_length_change=(len(s[a:b+1].replace('-',''))!=len(p['F']) or len(s[c:e+1].replace('-',''))!=len(p['R']))))
   seq=s[b+1:c].replace('-','')
   if dd[acc]['evaluable']=='1' and seq and set(seq)<=set('ACGT'):ins[acc]=seq
  if not name.endswith('D'):
   matrix={};ids=list(ins)
   for i,x in enumerate(ids):
    for y in ids[i+1:]:matrix[x,y]=matrix[y,x]=int(-aligner.score(ins[x],ins[y]))
   for x in ids:
    best=min(matrix[x,y] for y in ids if y!=x);nearest=[y for y in ids if y!=x and matrix[x,y]==best];identical=[y for y in ids if y==x or matrix[x,y]==0]
    genus=lambda z:meta[z]['organism'].split()[0]
    res['nearest_neighbors'].append(dict(gene=gene,window=name,accession=x,species=meta[x]['organism'],nearest_edit_distance=best,nearest_accessions=';'.join(nearest),genus_unique_at_identity=len({genus(y) for y in identical})==1,family_unique_at_identity=len({meta[y]['family'] for y in identical})==1,leave_one_genus_correct=all(genus(y)==genus(x) for y in nearest),leave_one_family_correct=all(meta[y]['family']==meta[x]['family'] for y in nearest)))
  for side in 'FR':
   variants=expand(p[side]);oligos.append((gene,name,side,variants))
   for scenario in ['50_per_variant','200_per_direction']:
    conc=50 if scenario=='50_per_variant' else 200/len(variants);tms=[primer3.calc_tm(v,dna_conc=conc,**kw) for v in variants];gc=[100*sum(c in 'GC' for c in v)/len(v) for v in variants]
    res['oligo_properties'].append(dict(gene=gene,pair=name,side=side,sequence=p[side],length=len(p[side]),variants=len(variants),scenario=scenario,variant_nM=conc,GC_min=min(gc),GC_max=max(gc),Tm_min=min(tms),Tm_max=max(tms)))
  if name.endswith('D'):
   op=next(q for q in ps if q['name']==name[:-1]);nested=all(set(I[x])<=set(I[y]) for side in 'FR' for x,y in zip(op[side],p[side]))
   if not nested:
    for fixed,deg in itertools.product([0,1],repeat=2):
     qs={side:''.join(y if (len(I[y])>1 and deg) or (len(I[y])==1 and fixed) else x for x,y in zip(op[side],p[side])) for side in 'FR'};cov=ex=0
     for acc in raw:
      if dd[acc]['evaluable']!='1':continue
      exact,covered,_=flags(raw[acc],qs['F'],qs['R'],True,False);ex+=exact;cov+=covered
     res['factorial_controls'].append(dict(gene=gene,pair=name,fixed_changes=fixed,degenerate_changes=deg,F=qs['F'],R=qs['R'],covered=cov,exact=ex,n=int(p['evaluable'])))
 print(gene,'extensions complete',flush=True)
for i,(g,n,s,vs) in enumerate(oligos):
 for h,m,t,ws in oligos[i+1:]:
  if g!=h or n==m:continue
  worst=max((primer3.calc_heterodimer(v,w,dna_conc=50,**kw).tm,v,w) for v in vs for w in ws)
  res['cross_assay_dimers'].append(dict(gene=g,pair1=n,end1=s,pair2=m,end2=t,max_Tm_50nM=worst[0],variant1=worst[1],variant2=worst[2],exceeds_55C=worst[0]>=55))
support={(r['gene'],r['pair'],r['side'],int(r['position'])):r for r in rows(ROOT/'audit/position_support.csv')};base={(r['gene'],r['name']):r for r in rows(ROOT/'all_results.csv')}
for t in rows(ROOT/'audit/consensus_trials.csv'):
 edits=re.findall(r'([FR])(\d+)([ACGT])>([ACGT])',t['edits']);freq=min(int(support[t['gene'],t['pair'],side,int(pos)][new])/int(support[t['gene'],t['pair'],side,int(pos)]['n']) for side,pos,old,new in edits);b=base[t['gene'],t['pair']]
 for dt in [4,5,6,7]:
  for fr in [.8,.85,.9,.95]:
   reasons=[]
   if float(t['max_delta_Tm'])>dt:reasons.append('Tm_difference')
   if freq<fr:reasons.append('frequency')
   if int(t['lost_covered']):reasons.append('lost_species')
   if int(t['added_background_model_hits']):reasons.append('new_background')
   if (int(t['exact'])-int(b['exact']))/int(b['evaluable'])<.1:reasons.append('exact_gain_below_10pp')
   if t['terminal_ok']!='True':reasons.append('terminal_filter')
   if float(t['max_structure_Tm'])>=55:reasons.append('structure')
   if float(t['min_Tm'])<45:reasons.append('minimum_Tm')
   res['threshold_sensitivity'].append(dict(gene=t['gene'],pair=t['pair'],edits=t['edits'],minimum_frequency=freq,delta_Tm_limit=dt,frequency_limit=fr,passes=not reasons,reasons=';'.join(reasons)))
for n,rr in res.items():write(n,rr)
print({n:len(rr) for n,rr in res.items()})
