from pathlib import Path
import sys,itertools
import numpy as np,pandas as pd
from scipy.stats import t

b=Path(__file__).resolve().parents[2]
d=pd.read_csv(b/'data/components/macrophages/analysis/coexpression_donor_contrasts.csv')
d=d[(d.cohort=='GSE248453')&(d.method=='depth_quartile_matched')].copy()
assert d.groupby('gene').size().eq(8).all() and d['gene'].nunique()==8
rows=[]
for gene,g in d.groupby('gene'):
 v=g.delta.to_numpy();n=len(v);avg=v.mean();se=v.std(ddof=1)/np.sqrt(n)
 signs=np.array(list(itertools.product([-1,1],repeat=n)))
 p=(np.abs(signs@v/n)>=abs(avg)-1e-12).mean()
 rows.append({'gene':gene,'n_donors':n,'mean':avg,'sd':v.std(ddof=1),'t_ci_low':avg-t.ppf(.975,n-1)*se,'t_ci_high':avg+t.ppf(.975,n-1)*se,'positive':sum(v>0),'negative':sum(v<0),'signflip_p':p,'min_two_sided_p':2/2**n})
out=pd.DataFrame(rows);order=np.argsort(out.signflip_p.to_numpy()); ps=out.signflip_p.to_numpy()[order]; adj=np.minimum.accumulate((ps*len(ps)/np.arange(1,len(ps)+1))[::-1])[::-1]; q=np.empty(len(ps));q[order]=np.minimum(adj,1);out['BH_q']=q
o=Path(sys.argv[1]);o.mkdir(parents=True,exist_ok=True)
d.to_csv(o/'Table3_GSE248453_逐供体四分位匹配.csv',index=False,encoding='utf-8-sig')
out.to_csv(o/'Table3_GSE248453_八基因供体汇总.csv',index=False,encoding='utf-8-sig')
print(out.round(4).to_string(index=False))

