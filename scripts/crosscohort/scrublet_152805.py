from pathlib import Path
import os,sys,gzip,traceback
os.environ['OPENBLAS_NUM_THREADS']='2';os.environ['OMP_NUM_THREADS']='2';os.environ['MKL_NUM_THREADS']='2'
import numpy as np,pandas as pd,scrublet as scr
from scipy.io import mmread
base=Path(__file__).resolve().parents[2];src=base/'data/external/GSE152805';out=Path(sys.argv[1])
md=pd.read_csv(src/'GSE152805_细胞聚类标注候选.csv.gz'); rows=[];summary=[]
for sample,gsm in [('SY_113','GSM4626763'),('SY_116','GSM4626764'),('SY_118','GSM4626765')]:
 x=md[md['sample'].eq(sample)];idx=x.index_zero_based.to_numpy();f=next((src/'滑膜原始矩阵').glob(gsm+'*.matrix.mtx.gz'))
 try:
  with gzip.open(f,'rb') as h: m=mmread(h).tocsc()
  m.sum_duplicates();raw=m[:,idx].T.tocsc();del m
  model=scr.Scrublet(raw,expected_doublet_rate=.06,sim_doublet_ratio=2,random_state=20261119)
  scores,pred=model.scrub_doublets(min_counts=2,min_cells=3,min_gene_variability_pctl=85,n_prin_comps=30,use_approx_neighbors=False,verbose=False)
  if pred is None:raise ValueError('automatic threshold not obtained')
  rows.extend([{'sample':sample,'index_zero_based':int(i),'doublet_score':float(s),'predicted_doublet':bool(b)} for i,s,b in zip(idx,scores,pred)])
  summary.append({'sample':sample,'total':len(idx),'predicted_doublets':int(sum(pred)),'fraction':float(np.mean(pred)),'automatic_threshold':float(model.threshold_),'expected_doublet_rate_assumed':.06,'status':'ok'})
 except Exception as e:traceback.print_exc();summary.append({'sample':sample,'status':'failed','error':repr(e)})
pd.DataFrame(rows).to_csv(out/'GSE152805_Scrublet_逐细胞.csv.gz',index=False,compression='gzip')
pd.DataFrame(summary).to_csv(out/'GSE152805_Scrublet_汇总.csv',index=False,encoding='utf-8-sig')
print(pd.DataFrame(summary).to_string(index=False))
if any(x['status']!='ok' for x in summary):sys.exit(1)

