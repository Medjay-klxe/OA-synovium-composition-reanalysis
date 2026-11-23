import sys
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from docx import Document
root=Path(__file__).resolve().parents[2]
base=root;out=base/'results/figures';out.mkdir(parents=True,exist_ok=True)
d=pd.read_csv(base/'data/umap'/'GSE216651_UMAP_7492.csv');lo=pd.read_csv(root/'data/components/macrophages/analysis'/'decomposition_LODO.csv').groupby('omitted_donor')[['composition','within_profile']].sum();assert len(d)==7492 and len(lo)==6
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],'font.size':7.2,'axes.labelsize':7.5,'xtick.labelsize':7,'ytick.labelsize':7,'pdf.fonttype':42,'svg.fonttype':'none'})
blue='#0072B2';orange='#D55E00';fig=plt.figure(figsize=(170/25.4,90/25.4),facecolor='white')
fig.text(.015,.98,'Fig. S1',fontsize=8,fontweight='bold',va='top')
a=fig.add_axes([.075,.31,.255,.51]);b=fig.add_axes([.46,.31,.20,.51]);c=fig.add_axes([.78,.31,.18,.51])
for grp,col in [('Healthy',blue),('OA',orange)]:
 q=d[d.group==grp];a.scatter(q.UMAP1,q.UMAP2,s=1.1,color=col,label=grp,alpha=.55,linewidths=0)
a.set_xlim(-12.5,7.5);a.set_ylim(-8,7);a.set_xticks([-10,-5,0,5]);a.set_yticks([-7.5,-5,-2.5,0,2.5,5]);a.set_aspect('equal');a.set_xlabel('UMAP 1');a.set_ylabel('UMAP 2');a.legend(loc='upper center',bbox_to_anchor=(.5,-.27),ncol=2,frameon=False,markerscale=4,fontsize=7)
b.barh(['Composition','Conditional profile'],[-.0246,.0992],color=[blue,orange]);b.invert_yaxis();b.axvline(0,color='black',lw=.6);b.set_xlabel('Two-component difference');b.set_xlim(-.05,.15);b.text(.008,0,'−0.0246',va='center',fontsize=7);b.text(.1042,1,'+0.0992',va='center',fontsize=7)
for i,k in enumerate(['composition','within_profile']): c.scatter(np.full(len(lo),i)+np.linspace(-.12,.12,len(lo)),lo[k],color=[blue,orange][i],s=18)
c.axhline(0,color='grey',lw=.6);c.set_xticks([0,1],['Composition','Conditional\nprofile']);c.set_ylabel('Two-component difference')
for x in [a,b,c]:x.spines[['top','right']].set_visible(False);x.tick_params(length=2)
for x,l in zip([.047,.43,.75],'ABC'):fig.text(x,.90,l,fontweight='bold',fontsize=11,va='center')
for x,t in zip([.2025,.56,.87],['Donor group','Two-component decomposition','Leave-one-donor-out']):fig.text(x,.90,t,ha='center',va='center',fontsize=7.5)
for ext in ['pdf','png','svg']:
 fig.savefig(out/('Supplementary_FigS1_供体组与分解敏感性.'+ext),dpi=300,facecolor='white')
 fig.savefig(out/('Additional_file_2.'+ext),dpi=300,facecolor='white')
plt.close(fig)