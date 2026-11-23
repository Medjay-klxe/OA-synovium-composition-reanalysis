from pathlib import Path
import os,sys
import pandas as pd,numpy as np
b=Path(__file__).resolve().parents[2];out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=True)
configs=[('GSE283080','data/external/GSE283080',[0,1,2,11,13]),('GSE152805','data/external/GSE152805',[8,13,17])]
rows=[]
for dataset,folder,clusters in configs:
 p=b/folder
 md=pd.read_csv(p/f'{dataset}_细胞聚类标注候选.csv.gz')
 z=np.load(p/f'{dataset}_聚类标记表达.npz')
 values=z['markers']; genes=list(z['genes']);assert len(md)==len(values)
 for sample in md['sample'].unique():
  whole=md['sample'].eq(sample).to_numpy();sel=whole & md.cluster.isin(clusters).to_numpy()
  subset=md.loc[sel]; v=values[sel];assert len(v)==len(subset)
  olr=v[:,genes.index('OLR1')]
  row={'dataset':dataset,'sample':sample,'qc_cells':int(whole.sum()),'candidate_cells':len(v),'candidate_share_among_qc':len(v)/whole.sum(),'umi_median':float(np.median(subset.total_umi)),'umi_q25':float(np.quantile(subset.total_umi,.25)),'umi_q75':float(np.quantile(subset.total_umi,.75))}
  for g in ['OLR1','GPNMB','MARCO','TREM2','APOE','LYVE1','FOLR2']:
   gene=v[:,genes.index(g)];row[g+'_detect_all']=np.mean(gene>0);row[g+'_mean_logcp10k']=np.mean(gene)
   for label,ix in [('OLR1_detected',olr>0),('OLR1_not_detected',olr==0)]:
    row[g+'_detect_'+label]=np.mean(gene[ix]>0) if ix.sum() else np.nan
  row['OLR1_positive_mean']=float(olr[olr>0].mean()) if (olr>0).any() else np.nan
  row['OLR1_positive_umi_median']=float(np.median(subset.total_umi.to_numpy()[olr>0])) if (olr>0).any() else np.nan
  row['OLR1_negative_umi_median']=float(np.median(subset.total_umi.to_numpy()[olr==0])) if (olr==0).any() else np.nan
  rows.append(row)
df=pd.DataFrame(rows);df.to_csv(out/'Fig6_候选髓系供体测序与检出指标.csv',index=False,encoding='utf-8-sig')
for dataset,folder,clusters in configs:
 e=pd.read_csv(b/folder/f'{dataset}_全转录组聚类匹配供体差值.csv')
 e['dataset']=dataset
 if dataset=='GSE283080': ef=e
 else: ef=pd.concat([ef,e],ignore_index=True)
ef.to_csv(out/'Fig6_跨队列相同方法供体效应.csv',index=False,encoding='utf-8-sig')
print(df[['dataset','sample','qc_cells','candidate_cells','candidate_share_among_qc','umi_median','OLR1_detect_all','GPNMB_detect_all','MARCO_detect_all']].round(3).to_string(index=False))
print(ef[ef.gene.isin(['GPNMB','MARCO'])][['dataset','sample','gene','effect']].round(3).to_string(index=False))
