from pathlib import Path
import sys,itertools
import numpy as np,pandas as pd
base=Path(__file__).resolve().parents[2];src=base/'data/external/GSE152805';out=Path(sys.argv[1])
md=pd.read_csv(src/'GSE152805_细胞聚类标注候选.csv.gz');z=np.load(src/'GSE152805_聚类标记表达.npz');X=z['markers'];genes=list(z['genes']);pred=pd.read_csv(out/'GSE152805_Scrublet_逐细胞.csv.gz')
assert len(md)==len(X) and len(pred)==len(md) and not pred.duplicated(['sample','index_zero_based']).any()
df=md.merge(pred[['sample','index_zero_based','predicted_doublet']],on=['sample','index_zero_based'],validate='one_to_one',sort=False);assert len(df)==len(md) and (df['sample'].to_numpy()==md['sample'].to_numpy()).all() and (df.index_zero_based.to_numpy()==md.index_zero_based.to_numpy()).all()
targets=['TREM2','GPNMB','APOE','MARCO','LYVE1','FOLR2','MRC1','F13A1'];idx=[genes.index(g) for g in targets];rng=np.random.default_rng(20261120);res=[];status=[]
for s in df['sample'].unique():
 mask=(df['sample'].eq(s)&df.cluster.isin([8,13,17])&~df.predicted_doublet).to_numpy();m=df.loc[mask].reset_index(drop=True);e=X[mask];pos=e[:,genes.index('OLR1')]>0
 edges=np.unique(np.quantile(m.total_umi,np.linspace(0,1,11)));bins=np.searchsorted(edges[1:-1],m.total_umi,side='right')+10*m.cluster.to_numpy();parts=[]
 for b in np.unique(bins):
  p=np.flatnonzero((bins==b)&pos);n=np.flatnonzero((bins==b)&~pos);k=min(len(p),len(n))
  if k:parts.append((p,n,k))
 n=sum(k for p,q,k in parts);row={'sample':s,'candidate_cells':len(m),'OLR1_detected':int(pos.sum()),'matched_pairs':n,'eligible':n>=30};status.append(row)
 if n<30:continue
 delta=[]
 for i in range(100):
  p=np.concatenate([rng.choice(a,k,False) for a,b,k in parts]);q=np.concatenate([rng.choice(b,k,False) for a,b,k in parts]);delta.append(e[p][:,idx].mean(axis=0)-e[q][:,idx].mean(axis=0))
 for j,g in enumerate(targets):res.append({'sample':s,'gene':g,'effect':np.mean(delta,axis=0)[j]})
pd.DataFrame(status).to_csv(out/'GSE152805_双细胞排除_匹配情况.csv',index=False,encoding='utf-8-sig')
r=pd.DataFrame(res);r.to_csv(out/'GSE152805_双细胞排除_供体差值.csv',index=False,encoding='utf-8-sig')
print(pd.DataFrame(status).to_string(index=False));print(r.groupby('gene').effect.agg(['mean',lambda v:sum(v>0)]).round(4).to_string())

