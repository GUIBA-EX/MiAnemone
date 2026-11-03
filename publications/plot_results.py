"""Generate manuscript figures from the archived assay and reference panels."""
import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import Counter

root=Path(__file__).resolve().parents[1]
def rows(p):
    with (root/p).open() as f:return list(csv.DictReader(f))
def save(fig,name):
    fig.tight_layout()
    fig.savefig(root/'publications/figures'/f'{name}.svg',bbox_inches='tight')
    plt.close(fig)
plt.rcParams.update({'svg.fonttype':'none','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
genes=['12S','16S']
fig,axs=plt.subplots(2,2,figsize=(9,5.8),sharey=True)
for i,g in enumerate(genes):
    rr={x['name']:x for x in rows(f'{g}/candidate_results.csv')}
    names=['XS','S','M','L','XL']+(['H'] if i else [])
    for ax,key in zip(axs[i],['covered','exact']):
        for j,(suffix,color) in enumerate([('','#277884'),('D','#c37325')]):
            ax.bar([k+(-.18 if not j else .18) for k in range(len(names))],[100*int(rr[n+suffix][key])/int(rr[n+suffix]['evaluable']) for n in names],.35,color=color,label=['Ordinary','D'][j])
        ax.set(ylim=(0,105),ylabel='% of evaluable species',title=g+' · '+key,xticks=range(len(names)),xticklabels=[f"{n}\nn={rr[n]['evaluable']}" for n in names])
axs[0,0].legend(fontsize=8)
save(fig,'paired_matches')
fig,axs=plt.subplots(1,2,figsize=(9,3.1))
for ax,g,file in zip(axs,genes,['barcode_identity','insert_resolution']):
    rr=rows(f'{g}/{file}.csv')
    vals=[(int(x.get('unique',x.get('identity_unique'))),int(x['unambiguous'])) for x in rr]
    bars=ax.bar([x['window'] for x in rr],[100*n/d for n,d in vals],color='#277884')
    ax.set(title=g+' · unique reference inserts',ylim=(0,72),ylabel='% unique');ax.bar_label(bars,labels=[f'{n}/{d}' for n,d in vals],fontsize=8)
save(fig,'insert_identity')
fig,axs=plt.subplots(1,2,figsize=(9,3.4))
for ax,g,pair in zip(axs,genes,['MD','LD']):
    def group(x):
        for name in ['Antipatharia','Scleractinia','Corallimorpharia','Zoantharia','Octocorallia','Hydrozoa']:
            if name in x['lineage'].split(';'):return name
        return 'Other'
    panel=Counter(group(r) for r in rows(f'{g}/background_panel.csv'))
    hits=rows(f'{g}/background_hits.csv')
    labs=list(sorted(panel))
    for i,(model,k,color) in enumerate([('WC','strict','#277884'),('GT','strict','#83b8ba'),('WC','relaxed','#c37325'),('GT','relaxed','#edc79f')]):
        v=Counter(group(r) for r in hits if r['name']==pair and r['model']==model and r[k]=='1')
        ax.barh([j+(i-1.5)*.19 for j in range(len(labs))],[v[n] for n in labs],.18,label=model+' '+k,color=color)
    ax.set(title=g+' '+pair,yticks=range(len(labs)),yticklabels=[f'{n} ({panel[n]})' for n in labs],xlabel='Background records matched')
axs[1].legend(fontsize=7)
save(fig,'background')
