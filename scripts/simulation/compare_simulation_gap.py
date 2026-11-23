from pathlib import Path
import sys
import pandas as pd,numpy as np
from scipy.stats import t
b=Path(__file__).resolve().parents[2]
f=b/'data/simulation'
x=pd.concat([pd.read_csv(f/f'fits_{i}.csv.gz') for i in range(1,5)])
x=x[(x.reference=='all6')&(x.scenario=='both')].copy()
real=pd.read_csv(b/'data/components/bulk/analysis'/'all_models_standardized.csv')
real=real[(real.cohort=='primary_QC44')&(real.model=='unadjusted')].copy()
assert len(real)==78
def residual_sd(beta,p):
 p=np.clip(np.asarray(p),np.finfo(float).tiny,1-1e-15)
 statistic=t.isf(p/2,42)
 return abs(beta)/statistic/np.sqrt(1/25+1/19)
real['resid_sd']=residual_sd(real.beta,real.p)
x['resid_sd']=residual_sd(x.beta_before,x.p_before)
z=x.groupby(['batch','rep','method']).agg(median_abs_beta=('beta_before',lambda v:np.median(abs(v))),median_residual_sd=('resid_sd','median')).reset_index()
actual=real.groupby('method').agg(median_abs_beta=('beta',lambda v:np.median(abs(v))),median_residual_sd=('resid_sd','median')).reset_index()
z.to_csv(Path(sys.argv[1])/'Fig3_模拟评分异质性_逐次.csv',index=False,encoding='utf-8-sig')
actual.to_csv(Path(sys.argv[1])/'Fig3_真实评分异质性.csv',index=False,encoding='utf-8-sig')
print('ACTUAL',actual.to_string(index=False));print('SIM',z.groupby('method')[['median_abs_beta','median_residual_sd']].median().to_string())

