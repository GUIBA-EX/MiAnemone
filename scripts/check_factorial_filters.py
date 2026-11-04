import csv,json
from pathlib import Path
from review_chemistry import terminal,thermal
r=Path(__file__).resolve().parents[1];out=[]
for x in csv.DictReader((r/'revision/factorial_controls.csv').open()):
 ok,details=thermal(x['F'],x['R']);minimum=min(min(y['F_tm_min'],y['R_tm_min']) for y in details);term=terminal(x['F']) and terminal(x['R']);delta=max(abs((y['F_tm_min']+y['F_tm_max']-y['R_tm_min']-y['R_tm_max'])/2) for y in details)
 out.append(dict(x,terminal_pass=term,thermal_pass=ok and minimum>=45,filter_pass=term and ok and minimum>=45,max_delta_Tm=delta,min_variant_Tm=minimum,max_structure_Tm=max(y['structure_max'] for y in details)))
with (r/'revision/factorial_filters.csv').open('w') as f:
 w=csv.DictWriter(f,out[0]);w.writeheader();w.writerows(out)
print(json.dumps([x for x in out if x['fixed_changes']=='0' and x['degenerate_changes']=='1'],indent=2))
