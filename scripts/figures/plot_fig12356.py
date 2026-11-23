import sys,os

from pathlib import Path
import numpy as np,pandas as pd
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
base=Path(__file__).resolve().parents[2]; out=base/'results/figures';out.mkdir(parents=True,exist_ok=True)
p=base/'data/components'; s=base/'data/gap';f=base/'data/crosscohort'
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],'font.size':7.2,'axes.labelsize':7.5,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'pdf.fonttype':42,'ps.fonttype':42,'axes.linewidth':.7,'savefig.dpi':300})
W=170/25.4;blue='#0072B2';orange='#D55E00';green='#009E73';grey='#777777';rng=np.random.default_rng(29)
def read(path):return pd.read_csv(path)
def label(ax,x='A'):ax.text(-.13,1.07,x,transform=ax.transAxes,weight='bold',fontsize=10,va='top',clip_on=False)
def save(fig,n):
 fig.savefig(out/f'Fig{n}.pdf',facecolor='white');fig.savefig(out/f'Fig{n}.png',dpi=300,facecolor='white');plt.close(fig)
def style(ax):
 ax.spines[['top','right']].set_visible(False);ax.tick_params(length=2.5);ax.grid(axis='x',alpha=.13)


plt.rcParams.update({'mathtext.fontset':'custom','mathtext.rm':'Arial','mathtext.it':'Arial:italic','mathtext.bf':'Arial:bold'})
# Fig 1: three-tier workflow with study roles; all wording from final manuscript.
fig,axs=plt.subplots(3,1,figsize=(W,5.0));fig.subplots_adjust(left=.065,right=.985,top=.934,bottom=.06,hspace=.27)
for ax in axs:ax.set_axis_off()
def box(ax,x,y,w,h,t,face):
 ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.012,rounding_size=.02',lw=.8,ec='#9aaab3',fc=face));ax.text(x+w/2,y+h/2,t,ha='center',va='center',fontsize=7.2)
for i,(title,items) in enumerate([('Tissue-level analysis',[('GSE89408\n25 non-OA + 19 OA\nBulk profiles','#eaf3f7'),('39 gene sets\nTwo scoring methods','#edf4ee'),('Semi-synthetic mixtures\n4 scenarios × 100','#f8f2e6')]),('Macrophage subsets',[('GSE216651\n3 healthy + 3 OA','#eaf3f7'),('7,492 macrophages\n8 subsets','#edf4ee'),('Donor-level decomposition\n3 components','#f8f2e6')]),('Within-donor expression',[('GSE248453\n8 OA; 7 eligible','#eaf3f7'),('GSE283080\n4 OA','#edf4ee'),('GSE152805\n3 OA','#f8f2e6')])]):
 ax=axs[i];ax.text(0,1.04,chr(65+i),fontsize=10,weight='bold',transform=ax.transAxes);ax.text(.055,1.04,title,fontsize=8,weight='bold',transform=ax.transAxes)
 for j,(t,c) in enumerate(items):box(ax,.06+j*.32,.15,.28,.67,t,c)
 if i==2:ax.text(.5,-.04,'Within-OA comparisons; no independent OA vs healthy validation',ha='center',va='top',transform=ax.transAxes,fontsize=7)
save(fig,1)

# Fig 2: compare the 39 matching set-level coefficients for 44 samples.
x=read(p/'bulk/analysis/all_models_standardized.csv');x=x[(x.cohort=='primary_QC44')&(x.model=='macrophage')]
fig,axs=plt.subplots(1,2,figsize=(W,3.35));fig.subplots_adjust(left=.10,right=.98,bottom=.20,top=.85,wspace=.33)
for i,(ax,method,title) in enumerate(zip(axs,['z_mean','ssgsea'],['z-mean','ssGSEA'])):
 q=x[x.method==method];assert len(q)==39
 sig=q.baseline_q<.05;ax.scatter(q.baseline_beta[~sig],q.beta[~sig],s=15,c=grey,alpha=.7);ax.scatter(q.baseline_beta[sig],q.beta[sig],s=18,c=orange)
 lim=[min(-.5,q.baseline_beta.min(),q.beta.min())-.12,max(q.baseline_beta.max(),q.beta.max())+.12];ax.plot(lim,lim,ls='--',lw=.7,color='black');ax.axhline(0,color='grey',lw=.5);ax.axvline(0,color='grey',lw=.5);ax.set(xlim=lim,ylim=lim,xlabel='Unadjusted OA coefficient',ylabel='Macrophage-adjusted coefficient');ax.set_title(title,fontsize=8,pad=8);ax.text(.04,.94,f"BH $\mathit{{q}}$ < 0.05: {sum(sig)} → {sum(q.q<.05)}\nMedian decrease: {np.median(q.delta_beta):.4f}",transform=ax.transAxes,va='top',fontsize=7);label(ax,chr(65+i));style(ax)
axs[1].legend(handles=[Line2D([],[],marker='o',ls='',color=orange,label='Unadjusted $\mathit{q}$ < 0.05'),Line2D([],[],marker='o',ls='',color=grey,label='Unadjusted $\mathit{q}$ ≥ 0.05')],loc='lower right',frameon=False)
save(fig,2)

# Fig 3: scenario summary; paired replicate-level medians; sensitivity.
r=read(s/'Fig3_配对每次模拟结果.csv');gap=read(s/'Fig3_跨度与背景敏感性汇总.csv')
fig=plt.figure(figsize=(W,5.25));gs=GridSpec(2,2,figure=fig,left=.10,right=.97,bottom=.10,top=.93,hspace=.56,wspace=.40);axs=[fig.add_subplot(gs[i,j]) for i in range(2) for j in range(2)]
ax=axs[0];scens=['no_change','composition_only','expression_only','both'];labs=['Neither','Composition','Expression','Both'];pos=np.arange(4)
for k,(meth,col) in enumerate([('z_mean',blue),('ssGSEA',orange)]):
 vals=[(r[(r.scenario==sc)&(r.method==meth)].n_sig_est==0).mean()*100 for sc in scens];ax.bar(pos+(k-.5)*.44,vals,.26,color=col,label='z-mean' if k==0 else 'ssGSEA')
 ax.text(pos[-1]+(k-.5)*.44,vals[-1]+2,f'{vals[-1]:.0f}/100',ha='center',fontsize=6.5)
ax.set_xticks(pos,labs,rotation=18);ax.set_ylabel('Replicates with zero sets (%)');ax.set_ylim(0,115);ax.legend(frameon=False,loc='upper left',bbox_to_anchor=(0,1.19),ncol=2);label(ax,'A');style(ax)
for j,(meth,col) in enumerate([('z_mean',blue),('ssGSEA',orange)],1):
 ax=axs[j];q=r[(r.scenario=='both')&(r.method==meth)];ax.scatter(q.median_beta_oracle,q.median_beta_est,s=12,color=col,alpha=.45,edgecolor='none');lo=min(q.median_beta_oracle.min(),q.median_beta_est.min());hi=max(q.median_beta_oracle.max(),q.median_beta_est.max());ax.plot([lo,hi],[lo,hi],'--',color='black',lw=.6);ax.set_xlabel('Known simulated RNA fraction\nmedian adjusted coefficient');ax.set_ylabel('Estimated fraction\nmedian adjusted coefficient');ax.set_title('z-mean' if meth=='z_mean' else 'ssGSEA',fontsize=8);ax.text(.04,.97,('Known vs est.: 0.1271 vs 0.1681' if meth=='z_mean' else 'Known vs est.: 0.1689 vs 0.2203'),transform=ax.transAxes,va='top',fontsize=7);(ax.text(-.16,1.10,'B',transform=ax.transAxes,weight='bold',fontsize=10,va='top',clip_on=False) if j==1 else label(ax,chr(65+j)));style(ax)
ax=axs[3]; scenarios=['both_baseline','both_half_contrast','both_high_contrast','both_donor_background'];names=['Baseline','Composition 0.5×','Composition 1.5×','Variable background'];q=gap[(gap.method=='ssGSEA')&(gap.reference=='all6')].set_index('scenario'); y=np.arange(4)
for i,k in enumerate(scenarios):
 z=q.loc[k];ax.errorbar(z.median_drop,i,xerr=[[z.median_drop-z.q025],[z.q975-z.median_drop]],fmt='o',color=blue,capsize=2,ms=4,lw=.8)
ax.axvline(.2804,ls='--',color=orange,lw=1,label='Observed 0.2804');ax.set_yticks(y,names);ax.invert_yaxis();ax.set_xlabel('Median coefficient decrease\n(central 95% simulation range)');ax.legend(frameon=False,loc='lower right');label(ax,'D');style(ax)
save(fig,3)

# Fig 5: measured preprocessing coefficients; factorial marginal effects; conditional AUC difference.
fig,axs=plt.subplots(1,3,figsize=(W,3.55));fig.subplots_adjust(left=.09,right=.98,bottom=.26,top=.84,wspace=.48)
q=read(p/'candidates/analysis/Fig3A_preprocessing_coefficients.csv');ax=axs[0]
prep=['raw_log2_counts','log2_CPM_plus1','TMM_logCPM'];print('preprocess values',q.preprocessing.unique())
for i,gene in enumerate(['OLR1','PYCARD']):
 for j,pre in enumerate(prep):
  sub=q[(q.gene==gene)&(q.preprocessing==pre)];
  if sub.empty:continue
  z=sub.iloc[0];xx=i*4+j;ax.plot([xx,xx],[z.beta,z.beta_adj],color=grey,lw=.8);ax.scatter(xx,z.beta,facecolors='white',edgecolors=blue if i==0 else orange,s=25,zorder=3);ax.scatter(xx,z.beta_adj,color=blue if i==0 else orange,s=22,zorder=3)
ax.set_xticks([0,1,2,4,5,6],['Raw','CPM','TMM']*2,rotation=45);ax.text(.19,1.03,'OLR1',transform=ax.transAxes,ha='center');ax.text(.8,1.03,'PYCARD',transform=ax.transAxes,ha='center');ax.axhline(0,color='grey',lw=.6);ax.set_ylabel('OA coefficient');ax.legend(handles=[Line2D([],[],marker='o',ls='',markerfacecolor='white',markeredgecolor='black',label='Unadjusted'),Line2D([],[],marker='o',ls='',color='black',label='Adjusted')],loc='upper left',bbox_to_anchor=(-.1,-.23),ncol=2,frameon=False,fontsize=7);label(ax,'A');style(ax)
ax=axs[1];q=read(p/'mixtures/analysis/factorial_primary_effects.csv').set_index('variable');xx=np.arange(2);ax.bar(xx-.17,[q.loc[g,'composition'] for g in ['OLR1','PYCARD']],.32,color=blue,label='Composition');ax.bar(xx+.17,[q.loc[g,'profile'] for g in ['OLR1','PYCARD']],.32,color=orange,label='Profile');ax.axhline(0,color='black',lw=.6);ax.set_xticks(xx,['OLR1','PYCARD']);ax.set_ylabel('Marginal contrast');ax.set_xlabel('7.6% vs 20.7%; healthy vs OA',fontsize=7);ax.legend(frameon=False,fontsize=7,loc='best');label(ax,'B');style(ax)
ax=axs[2];r=read(p/'candidates/analysis/CV_incremental_AUC_candidates.csv').set_index('gene').loc['OLR1'];ax.axvline(0,color='grey',ls='--',lw=.6);ax.errorbar(r.mean_delta_AUC,0,xerr=[[r.mean_delta_AUC-r.conditional_lo],[r.conditional_hi-r.mean_delta_AUC]],fmt='o',color=orange,capsize=3);ax.set_yticks([0],['+OLR1']);ax.set_xlabel('ΔAUC (conditional 95% CI)');ax.text(.03,.95,'AUC 0.799 → 0.816',transform=ax.transAxes,va='top');ax.set_xlim(-.065,.095);label(ax,'C');style(ax)
fig.text(.98,.025,'Fixed OOF predictions\nand composition estimates',ha='right',fontsize=7)
save(fig,5)

# Fig 6: 8 genes × 5 donor cohorts, donor means with CI. Sensitivity contrasts separate.
raw=read(f/'统一十等分_逐供体差值.csv');summ=read(f/'统一十等分_队列分组汇总.csv');quart=read(base/'data/crosscohort/Table3_GSE248453_八基因供体汇总.csv');dbl=read(f/'GSE152805_双细胞排除_供体差值.csv')
genes=['TREM2','GPNMB','APOE','MARCO','LYVE1','FOLR2','MRC1','F13A1']; groups=[('GSE216651','Healthy','GSE216651 healthy\n3/3'),('GSE216651','OA','GSE216651 OA\n3/3'),('GSE248453','OA','GSE248453 OA\n7/8'),('GSE283080','OA','GSE283080 OA\n4/4'),('GSE152805','OA','GSE152805 OA\n3/3')]
fig=plt.figure(figsize=(W,7.35));gs=GridSpec(2,5,figure=fig,height_ratios=[2.6,1],left=.13,right=.92,bottom=.15,top=.93,hspace=.43,wspace=.17)
for j,(ds,grp,title) in enumerate(groups):
 ax=fig.add_subplot(gs[0,j]);ax.axvline(0,color='grey',lw=.6);ax.set_title(title,fontsize=7);ax.set_xlim(-2.3,2.3);ax.set_ylim(7.6,-.6)
 for i,g in enumerate(genes):
  q=raw[(raw.dataset==ds)&(raw.group==grp)&(raw.gene==g)];z=summ[(summ.dataset==ds)&(summ.group==grp)&(summ.gene==g)].iloc[0]
  ax.scatter(q.effect,np.full(len(q),i)+np.linspace(-.15,.15,len(q)),s=9,color=grey,alpha=.75,zorder=2);ax.errorbar(z['mean'],i,xerr=[[z['mean']-z.ci_low],[z.ci_high-z['mean']]],fmt='o',ms=3.3,capsize=1.5,color=blue,lw=.8,zorder=4)
 ax.set_yticks(range(8),genes if j==0 else []);ax.set_xticks([-2,0,2]);ax.set_xlabel('Δ expression');ax.spines[['top','right']].set_visible(False)
 if j==0:ax.text(-.21,1.07,'A',transform=ax.transAxes,weight='bold',fontsize=10,va='top',clip_on=False)
ax=fig.add_subplot(gs[1,:3]); ax.set_title('GSE248453 OA',fontsize=7,pad=5); q=quart.set_index('gene').loc[genes];yy=np.arange(8);ax.axvline(0,color='grey',lw=.6);ax.errorbar(q['mean'],yy,xerr=[q['mean']-q.t_ci_low,q.t_ci_high-q['mean']],fmt='o',ms=3,color=green,capsize=2,lw=.8);ax.set_yticks(yy,genes);ax.invert_yaxis();ax.set_xlabel('Quartile matched, 8 OA donors');ax.text(.98,.03,'Six genes: BH $\mathit{q}$ = 0.0104',transform=ax.transAxes,ha='right',fontsize=7);label(ax,'B');style(ax)
ax=fig.add_subplot(gs[1,3:]);ax.set_title('GSE152805 OA',fontsize=7,pad=5);ax.axvline(0,color='grey',lw=.6)
for i,g in enumerate(['GPNMB','MARCO']):
 q=dbl[dbl.gene==g].effect;ax.scatter(q,np.full(len(q),i)+np.linspace(-.12,.12,len(q)),s=15,color=grey);ax.scatter(q.mean(),i,s=30,color=orange)
ax.set_yticks([0,1],['GPNMB','MARCO']);ax.yaxis.tick_right();ax.invert_yaxis();ax.set_xlabel('Cluster + UMI;\nScrublet excluded');ax.set_xlim(-.85,.1);label(ax,'C');style(ax)
fig.text(.13,.025,'Panel A: donors; mean and 95% CI',ha='left',fontsize=7);fig.text(.68,.025,'Gray: donors; orange: mean',ha='left',fontsize=7)
save(fig,6)
print('DONE',out)

