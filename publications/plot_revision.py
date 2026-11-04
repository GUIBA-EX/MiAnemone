import csv,math
from collections import defaultdict,Counter
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
font_manager.fontManager.addfont('/System/Library/Fonts/Avenir Next.ttc')
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'publications/figures'
def rows(p):return list(csv.DictReader((ROOT/p).open()))
def save(fig,n):
 for t in fig.findobj(matplotlib.text.Text): t.set_text(t.get_text().replace(chr(8242),chr(39)))
 fig.savefig(OUT/(n+'.svg'),bbox_inches='tight');fig.savefig(OUT/(n+'.png'),dpi=180,bbox_inches='tight');plt.close(fig)
def wilson(k,n):
 z=1.96;p=k/n;d=1+z*z/n;v=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n));return (p+z*z/(2*n)-v)/d,(p+z*z/(2*n)+v)/d
plt.rcParams.update({'font.family':'Avenir Next','font.size':10,'axes.titlesize':11,'axes.labelsize':9,'axes.linewidth':0.7,'xtick.labelsize':8,'ytick.labelsize':8,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
window_colors={'XS':'#159fa8','S':'#3d7993','M':'#ce962d','L':'#4779b8','XL':'#9063a8','H':'#82888c'}
colors=['#287c8e','#d58b34'];names=['XS','S','M','L','XL','H'];genes=['12S','16S']
fig,axs=plt.subplots(2,2,figsize=(9,6),layout='constrained')
for i,g in enumerate(genes):
 rr={r['pair']:r for r in rows('revision/panel_comparison.csv') if r['gene']==g and r['panel']=='common'};nn=names[:5]+(['H'] if i else [])
 for ax,key in zip(axs[i],['covered','exact']):
  for j,suffix in enumerate(['','D']):
   vals=[rr[n+suffix] for n in nn];ns=[int(v['n']) for v in vals];ks=[int(v[key]) for v in vals];ys=[100*k/n for k,n in zip(ks,ns)];ci=[wilson(k,n) for k,n in zip(ks,ns)];errs=[[y-lo*100 for y,(lo,hi) in zip(ys,ci)],[hi*100-y for y,(lo,hi) in zip(ys,ci)]]
   bars=ax.bar(np.arange(len(nn))+(j-.5)*.35,ys,.34,color=[('#82888c' if j==0 else '#b8bdc0') if n=='H' else colors[j] for n in nn],label=['Ordinary','D'][j]);ax.errorbar([b.get_x()+b.get_width()/2 for b in bars],ys,yerr=errs,fmt='none',ecolor='#333',elinewidth=.7,capsize=2)
   for b,k,y,limits in zip(bars,ks,ys,ci):ax.text(b.get_x()+b.get_width()/2,limits[1]*100+3,str(k),ha='center',fontsize=7)
  ax.set(xticks=range(len(nn)),xticklabels=nn,ylim=(0,117),yticks=[0,25,50,75,100],title=f'{g} · {key} · common n={ns[0]}',ylabel='Species (%)')
fig.legend(handles=[__import__('matplotlib.patches',fromlist=['Patch']).Patch(color=c,label=l) for c,l in zip(colors,['Ordinary','D'])],loc='outside upper center',ncol=2);save(fig,'review_common_panel')
# Separate entropy, called-species support and exact binding intervals.
fig=plt.figure(figsize=(9,8.4),layout='constrained')
subfigs=fig.subfigures(2,1,hspace=.06)
for sub,g in zip(subfigs,genes):
 axs=sub.subplots(3,1,sharex=True,gridspec_kw={'height_ratios':[1.25,.65,1.45]})
 rr=[r for r in rows('revision/entropy.csv') if r['gene']==g]
 x=np.array([int(r['coordinate']) for r in rr]);n=np.array([int(r['n']) for r in rr]);h=np.array([float(r['entropy']) if int(r['n'])>=2 else np.nan for r in rr])
 axs[0].plot(x,h,lw=.65,color='#527d70');axs[0].set(ylabel='Entropy (bits)',ylim=(0,2.05),yticks=[0,1,2],title=g)
 axs[1].fill_between(x,0,n,color='#cadcd5');axs[1].plot(x,n,lw=.65,color='#527d70');axs[1].set(ylabel='Called species',ylim=(0,280),yticks=[0,125,250])
 ps={p['name']:p for p in rows(g+'/candidate_results.csv')};nn=names[:5]+(['H'] if g=='16S' else [])
 for j,name in enumerate(nn):
  p=ps[name];d=ps[name+'D'];keys=['F_start','F_end','R_start','R_end'];assert all(p[k]==d[k] for k in keys)
  fs,fe,rs,re=[int(p[k]) for k in keys];col=window_colors[name]
  axs[2].plot([fe,rs],[j,j],lw=1,color=col,alpha=.5)
  # Arrow tip coincides with the inclusive 3-prime binding coordinate.
  axs[2].arrow(fs,j,fe-fs,0,width=.15,head_width=.34,head_length=min(7,fe-fs),length_includes_head=True,color=col)
  axs[2].arrow(re,j,rs-re,0,width=.15,head_width=.34,head_length=min(7,re-rs),length_includes_head=True,color=col)
 axs[2].set(yticks=range(len(nn)),yticklabels=[n+' / '+n+'D' for n in nn],ylim=(len(nn)-.5,-.7),xlabel='Gene-relative reference position (bp; 1-based)',xlim=(1,int(x[-1])))
 axs[2].text(.99,.99,'F: forward   R: reverse',transform=axs[2].transAxes,ha='right',va='top',fontsize=8)
 for ax in axs:ax.grid(axis='x',color='#dce4e1',lw=.4);ax.set_axisbelow(True)
save(fig,'review_entropy')
fig,axs=plt.subplots(2,2,figsize=(9,6),layout='constrained')
rr=rows('revision/unified_models.csv')
for i,g in enumerate(genes):
 for ax,model in zip(axs[i],['WC','GT']):
  pairs=[n+s for n in names[:5]+(['H'] if g=='16S' else []) for s in ['', 'D']];by={(r['pair'],r['role'],int(r['max_mismatches'])):r for r in rr if r['gene']==g and r['model']==model and r['terminal_zero']=='3'}
  for n in pairs:
   a=[by[n,'background',k] for k in range(5)];b=[by[n,'target',k] for k in range(5)];j=names.index(n.rstrip('D'));ax.plot([int(v['matched']) for v in a],[100*int(v['matched'])/int(v['n']) for v in b],marker='o',ms=2,lw=1,ls='--' if n.endswith('D') else '-',color=window_colors[n.rstrip('D')],label=n)
  ax.set(xlabel='Background records matched (count)',ylabel='Evaluable targets matched (%)',title=g+' · '+model,ylim=(0,105));ax.set_xscale('symlog',linthresh=1);ax.set_xlim(left=0);ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter());ax.legend(ncol=3,fontsize=7,loc='lower right')
save(fig,'review_unified_models')
fig,axs=plt.subplots(1,2,figsize=(9,4),layout='constrained')
for ax,g in zip(axs,genes):
 rr={r['pair']:r for r in rows('revision/joint_resolution.csv') if r['gene']==g and r['panel']=='common'};nn=names[:5]+(['H'] if g=='16S' else [])
 for j,suffix in enumerate(['','D']):
  v=[rr[n+suffix] for n in nn];ys=[100*int(r['covered_unique'])/int(r['n']) for r in v];bars=ax.bar(np.arange(len(nn))+(j-.5)*.35,ys,.34,color=[('#82888c' if j==0 else '#b8bdc0') if n=='H' else colors[j] for n in nn],label=['Ordinary','D'][j]);ax.bar_label(bars,labels=[r['covered_unique'] for r in v],fontsize=7,padding=2)
 ax.set(xticks=range(len(nn)),xticklabels=nn,ylabel='Covered and reference-unique (%)',title=g+f" · complete common inserts n={v[0]['n']}",ylim=(0,65));ax.legend(fontsize=8)
save(fig,'review_joint_resolution')
fig,axs=plt.subplots(2,2,figsize=(9,7),layout='constrained')
rs=rows('revision/mismatch_positions.csv')
for i,g in enumerate(genes):
 nn=[n+s for n in names[:5]+(['H'] if i else []) for s in ['','D']]
 for ax,side in zip(axs[i],['F','R']):
  data=defaultdict(list)
  for r in rs:
   if r['gene']==g and r['side']==side and r['gap']=='False' and r['template'] in 'ACGT':data[r['pair'],int(r['distance_from_3prime'])].append(int(r['mismatch']))
  mat=np.array([[100*np.mean(data[n,j]) if data[n,j] else np.nan for j in range(22)] for n in nn]);im=ax.imshow(mat,aspect='auto',vmin=0,vmax=100,cmap='magma_r');ax.set(yticks=range(len(nn)),yticklabels=nn,xticks=[0,2,4,9,14,19,21],title=g+' · '+side,xlabel="Distance from 3′ terminal base (0 = terminal)");ax.axvline(2.5,color='#167b93',lw=1)
fig.colorbar(im,ax=axs,shrink=.6,label='Mismatch among called projected sites (%)');save(fig,'review_mismatch_profiles')
# All original background outcomes, grouped consistently by accession lineage.
fig,axs=plt.subplots(2,1,figsize=(10,9),layout='constrained')
for ax,g in zip(axs,genes):
 meta={r['accession']:r for r in rows(g+'/background_panel.csv')};groups=['Scleractinia','Antipatharia','Corallimorpharia','Zoantharia','Octocorallia','Hydrozoa'];nn=[n+s for n in names[:5]+(['H'] if g=='16S' else []) for s in ['','D']];mat=np.zeros((len(nn),24),int)
 for r in rows(g+'/background_hits.csv'):
  tax=meta[r['accession']].get('lineage','');group=next((x for x in groups if x in tax.split(';')),None)
  if group is None:continue
  j=groups.index(group);m=0 if r['model']=='WC' else 1
  mat[nn.index(r['name']),j*4+m]+=int(r['strict']);mat[nn.index(r['name']),j*4+m+2]+=int(r['relaxed'])
 im=ax.imshow(mat,aspect='auto',cmap='Blues',norm=matplotlib.colors.SymLogNorm(linthresh=1,vmin=0,vmax=1000));ax.set(yticks=range(len(nn)),yticklabels=nn,xticks=np.arange(6)*4+1.5,xticklabels=groups,title=g);ax.tick_params(axis='x',labelrotation=60)
 for i in range(len(nn)):
  for j in range(24):
   if mat[i,j]:ax.text(j,i,str(mat[i,j]),ha='center',va='center',fontsize=8,color='white' if mat[i,j]>30 else 'black')
 ax.set_xticks(np.arange(24),minor=True);ax.set_xticklabels(['W','G','W+','G+']*6,minor=True,fontsize=7);ax.tick_params(axis='x',which='major',pad=17);ax.tick_params(axis='x',which='minor',length=0);
 for boundary in np.arange(3.5,23.5,4):ax.axvline(boundary,color='#b7c8cf',lw=.6)
 ax.set_xlabel('W / G: strict WC / GT; W+ / G+: relaxed WC / GT')
fig.colorbar(im,ax=axs,shrink=.65,pad=.02,ticks=[0,1,10,100,1000],label='Background records (symlog scale)');save(fig,'review_background_all')
