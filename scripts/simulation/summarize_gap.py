import sys
from pathlib import Path
import pandas as pd,numpy as np
p=Path(sys.argv[1]); d=pd.concat([pd.read_csv(p/f'fits_{i}.csv.gz') for i in [11,12,13,14]],ignore_index=True)
assert len(d)==100*4*2*2*39
d['drop']=d.beta_before-d.beta_after
r=d.groupby(['batch','rep','scenario','method','reference']).agg(drop=('drop','median'),sig=('q_after',lambda x:(x<.05).sum())).reset_index()
s=r.groupby(['scenario','method','reference']).agg(n=('rep','size'),median_drop=('drop','median'),q025=('drop',lambda x:x.quantile(.025)),q975=('drop',lambda x:x.quantile(.975)),mean_sig=('sig','mean'),zero_pct=('sig',lambda x:100*(x==0).mean())).reset_index()
s.to_csv(p/'Fig3_跨度与背景敏感性汇总.csv',index=False,encoding='utf-8-sig')
r.to_csv(p/'Fig3_跨度与背景_逐次汇总.csv',index=False,encoding='utf-8-sig')
print(s[(s.reference=='all6')&(s.method=='ssGSEA')].round(4).to_string(index=False))

