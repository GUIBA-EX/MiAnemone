import csv,collections,statistics
from pathlib import Path
from verify import rows,fasta
ROOT=Path(__file__).resolve().parents[1];summary=[];details=[];lengths=[]
for gene in ['12S','16S']:
 d=ROOT/'revision'/('all_records_'+gene);meta={r['accession']:r for r in rows(d/'panel.csv')};al={k.removeprefix('_R_'):s for k,s in fasta(d/'aligned.fasta').items()};co=[i for i,c in enumerate(al['NC_000933.1']) if c!='-'];hits={(r['name'],r['accession']):r for r in csv.DictReader((d/'independent_hits.tsv').open(),delimiter='\t') if r['role']=='target' and r['model']=='WC'}
 for p in rows(ROOT/gene/'candidate_results.csv'):
  a,b,c,e=[co[int(p[k])-1] for k in ['F_start','F_end','R_start','R_end']];groups=collections.defaultdict(list)
  for k,r in meta.items():
   s=al[k];ok=bool(s[:a+1].strip('-') and s[e:].strip('-') and set(s[a:b+1]+s[c:e+1])<=set('ACGT-'))
   if ok:groups[r['organism']].append((k,hits.get((p['name'],k),{}),s[a:b+1].replace('-','')+'|'+s[c:e+1].replace('-','')))
  for species,rs in groups.items():
   n=len(rs);cov=sum(int(r.get('le2last3',0)) for k,r,s in rs);ex=sum(int(r.get('exact',0)) for k,r,s in rs);details.append(dict(gene=gene,pair=p['name'],species=species,evaluable_records=n,covered_records=cov,exact_records=ex,distinct_binding_haplotypes=len({s for k,r,s in rs}),accessions=';'.join(k for k,r,s in rs)))
  rr=[r for r in details if r['gene']==gene and r['pair']==p['name']];summary.append(dict(gene=gene,pair=p['name'],species=len(rr),records=sum(r['evaluable_records'] for r in rr),any_covered=sum(r['covered_records']>0 for r in rr),majority_covered=sum(r['covered_records']>r['evaluable_records']/2 for r in rr),all_covered=sum(r['covered_records']==r['evaluable_records'] for r in rr),any_exact=sum(r['exact_records']>0 for r in rr),all_exact=sum(r['exact_records']==r['evaluable_records'] for r in rr),multirecord_species=sum(r['evaluable_records']>1 for r in rr),variable_binding_species=sum(r['distinct_binding_haplotypes']>1 for r in rr)))
for r in rows(ROOT/'all_results.csv'):
 rr=[int(x['product_bp']) for x in rows(ROOT/'revision/product_lengths.csv') if x['gene']==r['gene'] and x['pair']==r['name']];lengths.append(dict(gene=r['gene'],pair=r['name'],products=len(rr),min=min(rr),median=statistics.median(rr),max=max(rr),overlap20_2x150=sum(n<=280 for n in rr),overlap20_2x250=sum(n<=480 for n in rr),overlap20_2x300=sum(n<=580 for n in rr)))
for name,rr in [('all_record_summary',summary),('all_record_species',details),('length_summary',lengths)]:
 with (ROOT/'revision'/(name+'.csv')).open('w') as h:
  w=csv.DictWriter(h,rr[0]);w.writeheader();w.writerows(rr)
