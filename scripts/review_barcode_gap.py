"""Observed within-species distances and leave-one-record nearest-neighbour diagnostics."""
from collections import defaultdict
from Bio.Align import PairwiseAligner
from verify import rows,fasta
from pathlib import Path
import csv
ROOT=Path(__file__).resolve().parents[1];out=[];aln=PairwiseAligner(mode='global',match_score=0,mismatch_score=-1,open_gap_score=-1,extend_gap_score=-1)
for gene in ['12S','16S']:
 d=ROOT/'revision'/('all_records_'+gene);meta={r['accession']:r for r in rows(d/'panel.csv')};al={k.removeprefix('_R_'):v for k,v in fasta(d/'aligned.fasta').items()};co=[i for i,c in enumerate(al['NC_000933.1']) if c!='-']
 for p in rows(ROOT/gene/'candidate_results.csv'):
  if p['name'].endswith('D'):continue
  a,b=co[int(p['F_end'])-1],co[int(p['R_start'])-1];seq={};groups=defaultdict(list)
  for k in meta:
   s=al[k][a+1:b].replace('-','')
   if al[k][:a+1].strip('-') and al[k][b:].strip('-') and s and set(s)<=set('ACGT'):seq[k]=s;groups[meta[k]['organism']].append(k)
  cache={}
  def dist(k,j):
   x,y=sorted((seq[k],seq[j]));key=x,y
   if key not in cache:cache[key]=int(-aln.score(x,y))
   return cache[key]
  for sp,ids in groups.items():
   other=[k for k in seq if meta[k]['organism']!=sp];within=[dist(k,j) for i,k in enumerate(ids) for j in ids[i+1:]];between=min(dist(k,j) for k in ids for j in other)
   correct=ambiguous=0
   for k in ids:
    best=min(dist(k,j) for j in seq if j!=k);labels={meta[j]['organism'] for j in seq if j!=k and dist(k,j)==best};correct+=labels=={sp};ambiguous+=len(labels)>1
   out.append(dict(gene=gene,window=p['name'],species=sp,n=len(ids),within_max=max(within) if within else '',nearest_other_min=between,observed_gap=(max(within)<between) if within else '',leave_one_species_correct=correct,leave_one_ambiguous=ambiguous))
  print(gene,p['name'],len(seq),flush=True)
with (ROOT/'revision/observed_barcode_gap.csv').open('w') as h:
 w=csv.DictWriter(h,out[0]);w.writeheader();w.writerows(out)
