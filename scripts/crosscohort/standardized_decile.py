from pathlib import Path
import os,sys,itertools
os.environ['OPENBLAS_NUM_THREADS']='2'
import numpy as np,pandas as pd
from scipy.stats import t
from math import factorial
base=Path(__file__).resolve().parents[2];out=Path(sys.argv[1]);genes=['TREM2','GPNMB','APOE','MARCO','LYVE1','FOLR2','MRC1','F13A1'];rng=np.random.default_rng(20261126)
records=[];status=[]
for ds in ['GSE216651','GSE248453','GSE283080','GSE152805']:
 if ds in ['GSE216651','GSE248453']:
  md=pd.read_csv(out/f'{ds}_metadata.csv.gz');v=pd.read_csv(out/f'{ds}_logCP10K.csv.gz');assert np.array_equal(md.row_index,v.row_index)
  md=md.rename(columns={'library_total':'total_umi'});v=v[['OLR1']+genes].to_numpy()
 else:
  folder=base/('data/external/GSE283080' if ds=='GSE283080' else 'data/external/GSE152805')
  md=pd.read_csv(folder/f'{ds}_细胞聚类标注候选.csv.gz');z=np.load(folder/f'{ds}_聚类标记表达.npz');gs=list(z['genes']);v=z['markers'][:,[gs.index(g) for g in ['OLR1']+genes]]
  keep=md.cluster.isin([0,1,2,11,13] if ds=='GSE283080' else [8,13,17]).to_numpy();md=md.loc[keep].reset_index(drop=True);v=v[keep]
 for sample in md['sample'].unique():
  keep=md['sample'].eq(sample).to_numpy(); m=md.loc[keep].reset_index(drop=True);x=v[keep];positive=x[:,0]>0
  edges=np.unique(np.quantile(m.total_umi,np.linspace(0,1,11))); bins=np.searchsorted(edges[1:-1],m.total_umi,side='right')
  parts=[]
  for z in np.unique(bins):
   a=np.flatnonzero((bins==z)&positive);b=np.flatnonzero((bins==z)&~positive);n=min(len(a),len(b))
   if n:parts.append((a,b,n))
  pairs=sum(n for a,b,n in parts);group=('Healthy' if 'healthy' in sample else 'OA')
  row={'dataset':ds,'group':group,'sample':sample,'candidate_cells':len(m),'OLR1_detected':int(positive.sum()),'OLR1_undetected':int((~positive).sum()),'matched_pairs':pairs,'eligible':pairs>=30}
  if pairs>=30:
   diff=[];dep=[]
   for k in range(100):
    a=np.concatenate([rng.choice(a,n,False) for a,b,n in parts]);b=np.concatenate([rng.choice(b,n,False) for a,b,n in parts]);diff.append(x[a,1:].mean(axis=0)-x[b,1:].mean(axis=0));dep.append(np.log1p(m.total_umi.to_numpy()[a]).mean()-np.log1p(m.total_umi.to_numpy()[b]).mean())
   row['matched_logUMI_difference']=np.mean(dep);vals=np.array(diff)
   for j,g in enumerate(genes):records.append({'dataset':ds,'group':group,'sample':sample,'gene':g,'effect':vals[:,j].mean(),'repetition_sd':vals[:,j].std(ddof=1),'matched_pairs':pairs})
  status.append(row)
a=pd.DataFrame(records);q=[]
for (ds,group,gene),g in a.groupby(['dataset','group','gene']):
 v=g.effect.to_numpy();n=len(v);mean=v.mean();se=v.std(ddof=1)/np.sqrt(n) if n>1 else np.nan
 signs=np.array(list(itertools.product([-1,1],repeat=n)));p=np.mean(np.abs(signs@v/n)>=abs(mean)-1e-12)
 q.append({'dataset':ds,'group':group,'gene':gene,'n':n,'mean':mean,'ci_low':mean-t.ppf(.975,n-1)*se if n>1 else np.nan,'ci_high':mean+t.ppf(.975,n-1)*se if n>1 else np.nan,'positive':int((v>0).sum()),'negative':int((v<0).sum()),'p_signflip':p,'min_two_sided_p':2/(2**n)})
q=pd.DataFrame(q)
def bh(p):
 a=p.to_numpy();i=np.argsort(a);v=np.minimum.accumulate((a[i]*len(a)/np.arange(1,len(a)+1))[::-1])[::-1];z=np.empty(len(a));z[i]=np.minimum(1,v);return z
q['q_BH']=q.groupby(['dataset','group']).p_signflip.transform(bh)
for frame,name in [(pd.DataFrame(status),'匹配纳入情况.csv'),(a,'统一十等分_逐供体差值.csv'),(q,'统一十等分_队列分组汇总.csv')]:frame.to_csv(out/name,index=False,encoding='utf-8-sig')
print(pd.DataFrame(status).to_string(index=False));print(q.round(4).to_string(index=False))

