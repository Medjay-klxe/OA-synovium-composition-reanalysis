"""Run a defined stage from the released data and scripts."""
from pathlib import Path
import argparse,subprocess,sys,os,shutil
ROOT=Path(__file__).resolve().parent
os.environ['OA_REPO_ROOT']=str(ROOT)
def run(*args):
    subprocess.run([str(x) for x in args],cwd=ROOT,check=True)
def python(script,*args):run(sys.executable,ROOT/'scripts'/script,*args)
def r(script,*args):run(os.environ.get('RSCRIPT','Rscript'),ROOT/'scripts'/script,*args)
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=['figures','matching','gse152805','bulk','candidates','mixtures','decomposition','simulation','gap','external_bulk'])
    a=parser.parse_args();out=ROOT/'results'/a.stage;out.mkdir(parents=True,exist_ok=True)
    if a.stage=='figures':
        for s in ['figures/plot_fig12356.py','figures/plot_fig4.py','figures/plot_s1.py']:python(s)
    elif a.stage=='matching':
        for f in (ROOT/'data/crosscohort').glob('GSE*_*.csv.gz'):shutil.copy2(f,out/f.name)
        python('crosscohort/standardized_decile.py',out)
    elif a.stage=='gse152805':
        raw=ROOT/'data/external/GSE152805'/'滑膜原始矩阵'
        if not raw.exists():parser.error('Download the three GSE152805 synovial matrices into data/external/GSE152805/滑膜原始矩阵; see README.')
        python('crosscohort/scrublet_152805.py',out);python('crosscohort/filter_152805_doublets.py',out)
    elif a.stage=='external_bulk':r('analysis/external_bulk.R',out)
    elif a.stage=='bulk':r('analysis/bulk_gene_sets.R',ROOT/'data/components/bulk/inputs',out)
    elif a.stage=='candidates':r('analysis/candidate_models.R',ROOT/'data/components/candidates/inputs',out)
    elif a.stage=='mixtures':r('analysis/controlled_mixtures.R',ROOT/'data/components/mixtures/inputs',out)
    elif a.stage=='decomposition':
        r('analysis/macrophage_scores_and_two_component.R',ROOT/'data/components/macrophages/inputs',out)
        r('analysis/three_component_decomposition.R')
    elif a.stage in ('simulation','gap'):
        script='paired_scores.R' if a.stage=='simulation' else 'gap_sensitivity.R'
        batches=range(1,5) if a.stage=='simulation' else range(11,15)
        for batch in batches:r('simulation/'+script,ROOT/'scripts/simulation/simulation_setup.R',out,batch,25)
        if a.stage=='simulation':r('simulation/summarize.R',out)
        else:python('simulation/summarize_gap.py',out)
