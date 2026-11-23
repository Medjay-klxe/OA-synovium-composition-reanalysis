from pathlib import Path
import sys
import pandas as pd,numpy as np
root=Path(__file__).resolve().parents[2]
p=root/'data/simulation'
out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=True)
x=pd.concat([pd.read_csv(p/f'fits_{i}.csv.gz') for i in range(1,5)],ignore_index=True)
keys=['batch','rep','scenario','method','gene_set']
a=x[x.reference=='all6'];o=x[x.reference=='oracle']
m=a.merge(o[keys+['beta_after','q_after']],on=keys,validate='one_to_one',suffixes=('_est','_oracle'))
assert len(m)==400*2*39
m['estimated_minus_oracle']=m.beta_after_est-m.beta_after_oracle
m['abs_difference']=abs(m.estimated_minus_oracle)
m['abs_effect_change']=abs(m.beta_after_est)-abs(m.beta_after_oracle)
r=m.groupby(['batch','rep','scenario','method']).agg(median_est_minus_oracle=('estimated_minus_oracle','median'),median_abs_difference=('abs_difference','median'),median_abs_effect_change=('abs_effect_change','median'),n_sig_est=('q_after_est',lambda s:(s<.05).sum()),n_sig_oracle=('q_after_oracle',lambda s:(s<.05).sum()),median_beta_est=('beta_after_est','median'),median_beta_oracle=('beta_after_oracle','median')).reset_index()
summary=r.groupby(['scenario','method']).agg(runs=('rep','size'),median_signed_error=('median_est_minus_oracle','median'),median_abs_error=('median_abs_difference','median'),median_abs_change=('median_abs_effect_change','median'),median_est=('median_beta_est','median'),median_oracle=('median_beta_oracle','median'),mean_sig_est=('n_sig_est','mean'),mean_sig_oracle=('n_sig_oracle','mean')).reset_index()
m.to_csv(out/'Fig3_估计组成与已知组成_逐集合.csv.gz',index=False,compression='gzip')
r.to_csv(out/'Fig3_配对每次模拟结果.csv',index=False,encoding='utf-8-sig')
summary.to_csv(out/'Fig3_模拟效应对照汇总.csv',index=False,encoding='utf-8-sig')
print(summary.to_string(index=False))
data=pd.read_csv(root/'data/components/bulk/analysis'/'sample_inclusion_and_metadata.csv')
print('COLUMNS',list(data.columns));print(data.head(2).to_string(index=False))

