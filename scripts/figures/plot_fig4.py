import sys
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from docx import Document
root=Path(__file__).resolve().parents[2]
base=root
out=base/'results/figures';out.mkdir(parents=True,exist_ok=True)
d=pd.read_csv(base/'data/umap'/'GSE216651_UMAP_7492.csv');p=root/'data/components/macrophages/analysis'
don=pd.read_csv(p/'GSE216651_donor_scores.csv');lo=pd.read_csv(p/'decomposition_LODO.csv').groupby('omitted_donor')[['composition','within_profile']].sum()
assert len(d)==7492 and len(don)==6 and d.cell.is_unique
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],'font.size':7.4,'axes.labelsize':7.5,'xtick.labelsize':7,'ytick.labelsize':7,'pdf.fonttype':42,'svg.fonttype':'none'})
plt.rcParams.update({'mathtext.fontset':'custom','mathtext.rm':'Arial','mathtext.it':'Arial:italic','mathtext.bf':'Arial:bold'})
blue='#0072B2';orange='#D55E00';grey='#777777';green='#009E73';colors=['#0072B2','#D55E00','#009E73','#CC79A7','#E69F00','#56B4E9','#6D5A91','#777777'];levels=sorted(d.subtype.unique());cmap=LinearSegmentedColormap.from_list('olr1_gray',['#dedede','#563171','#b73779','#f16b50','#fdbb64','#fcfdbf'],N=256);lim=max(float(d.OLR1.quantile(.99)),1e-6)
def tag(ax,c):ax.text(-.16,1.05,c,transform=ax.transAxes,fontweight='bold',fontsize=10,va='bottom')
def style(ax):ax.spines[['top','right']].set_visible(False);ax.tick_params(length=2)
def umap(ax):
 ax.set_xlim(-12.5,7.5);ax.set_ylim(-8,7);ax.set_xticks([-10,-5,0,5]);ax.set_yticks([-7.5,-5,-2.5,0,2.5,5]);ax.set_aspect('equal');ax.set_xlabel('UMAP 1');ax.set_ylabel('UMAP 2');style(ax)
def save(f,n):
 for ext in ['pdf','png','svg']:f.savefig(out/('Fig4.'+ext),dpi=300,facecolor='white')
 plt.close(f)
fig=plt.figure(figsize=(170/25.4,155/25.4));a=fig.add_axes([.115,.66,.32,.26]);b=fig.add_axes([.58,.66,.32,.26]);c=fig.add_axes([.12,.18,.30,.25]);e=fig.add_axes([.63,.18,.32,.25])
for i,k in enumerate(levels):
 q=d[d.subtype==k];a.scatter(q.UMAP1,q.UMAP2,color=colors[i],s=1.1,label=k,linewidths=0)
umap(a);tag(a,'A');a.set_title('Macrophage subsets',fontsize=8);a.legend(loc='upper center',bbox_to_anchor=(.5,-.25),ncol=2,frameon=False,fontsize=7,markerscale=3,columnspacing=.6,handletextpad=.3)
z=d.sort_values('OLR1');im=b.scatter(z.UMAP1,z.UMAP2,c=z.OLR1,s=1.1,cmap=cmap,vmin=0,vmax=lim,linewidths=0);umap(b);tag(b,'B');b.set_title('OLR1 expression',fontsize=8)
pos=b.get_position();cax=fig.add_axes([.92,pos.y0,.012,pos.height]);cb=fig.colorbar(im,cax=cax);cb.solids.set_rasterized(False);cb.ax.tick_params(labelsize=7,length=2);cb.ax.set_title('OLR1\n(log-norm)',fontsize=7)
for i,(g,col) in enumerate([('Healthy',blue),('OA',orange)]):
 v=don.loc[don.group==g,'Composite'];c.scatter(np.full(len(v),i)+np.linspace(-.05,.05,len(v)),v,color=col,s=23);c.plot([i-.14,i+.14],[v.mean()]*2,color='black',lw=1)
c.set_xticks([0,1],['Healthy','OA']);c.set_ylabel('Donor composite score');c.set_title('Donor-level score',fontsize=8);tag(c,'C');style(c)
vals=[-.0178,.0974,-.0050,.0746];pv=[.20,.20,.50,.40];names=['Composition','Expression','Covariance','Total'];e.barh(names,vals,color=[blue,orange,grey,green]);e.invert_yaxis();e.axvline(0,color='black',lw=.6);e.set_xlim(-.04,.20);e.set_xlabel('Three-component difference');e.set_title('Decomposition',fontsize=8)
for i,(v,pv0) in enumerate(zip(vals,pv)):e.text(.105,i,f'{v:+.4f}; $\mathit{{P}}$ = {pv0:.2f}',va='center',fontsize=7)
tag(e,'D');style(e);save(fig,'Fig4_巨噬细胞UMAP与三项分解')