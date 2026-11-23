# Fig2 revision: explicit QC, TMM logCPM, matched-sample comparisons.
# Run from the delivery root: Rscript code/fig2_revision_analysis.R inputs analysis
# Outputs never overwrite the original project. Gene sets are fixed in inputs/.
args <- commandArgs(trailingOnly=TRUE)
inp <- if(length(args)>=1)args[1] else 'inputs'
out <- if(length(args)>=2)args[2] else 'analysis'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
suppressPackageStartupMessages({library(data.table);library(edgeR);library(GSVA)})
lines <- readLines(gzfile(file.path(inp,'GSE89408_counts.txt.gz')),warn=FALSE)
cnt <- fread(text=paste(lines,collapse='\n'));rm(lines)
genes <- cnt[[1]];M <- as.matrix(cnt[,-1]);storage.mode(M)<-'double';rownames(M)<-genes
# Count aggregation must sum duplicate gene identifiers (none expected here).
duplicate_count <- sum(duplicated(genes))
if(duplicate_count>0)M<-rowsum(M,group=rownames(M),reorder=FALSE)
stopifnot(all(is.finite(M)),all(M>=0),!anyDuplicated(colnames(M)))
meta<-as.data.frame(fread(file.path(inp,'deconv_proportions_218.csv')))
names(meta)[names(meta)=='sample']<-'sample_id'
meta<-meta[,c('sample_id','disease','age','sex','fit_r','Macrophage','Fibroblast')]
stopifnot(!anyDuplicated(meta$sample_id),all(colnames(M)%in%meta$sample_id))
meta<-meta[match(colnames(M),meta$sample_id),]
meta$OA<-as.integer(meta$disease=='Osteoarthritis')
meta$sex<-factor(meta$sex)
eligible<-meta$disease%in%c('Normal','Osteoarthritis')
meta50<-meta[eligible,];M50<-M[,eligible,drop=FALSE]
qc<-is.finite(meta50$fit_r)&meta50$fit_r>=.6
meta44<-meta50[qc,];M44<-M50[,qc,drop=FALSE]
stopifnot(nrow(meta50)==50,nrow(meta44)==44,sum(meta44$OA)==19)
ll<-strsplit(readLines(file.path(inp,'macrophage_39_c5.go.v2022.1.Hs.symbols.gmt'),warn=FALSE),'\t')
sets<-setNames(lapply(ll,function(x)unique(x[-c(1,2)])),vapply(ll,`[`,character(1),1))
stopifnot(length(sets)==39)
names4<-c('GOBP_MACROPHAGE_ACTIVATION','GOBP_MACROPHAGE_COLONY_STIMULATING_FACTOR_PRODUCTION','GOBP_MACROPHAGE_FUSION','GOBP_MACROPHAGE_INFLAMMATORY_PROTEIN_1_ALPHA_PRODUCTION')
coverage<-data.table(geneset=names(sets),n_defined=lengths(sets),n_in_counts=vapply(sets,function(x)sum(x%in%rownames(M)),integer(1)))
# Primary TMM factors fitted on the 44 retained OA/normal samples.
d44<-DGEList(M44)
keep<-filterByExpr(d44,group=meta44$disease)
d44<-d44[keep,,keep.lib.sizes=FALSE]
d44<-normLibSizes(d44,method='TMM')
expr44<-cpm(d44,log=TRUE,prior.count=1)
used<-lapply(sets,intersect,y=rownames(expr44))
coverage[,n_scored:=lengths(used)]
coverage[,score_eligible:=n_scored>=3]
fwrite(coverage,file.path(out,'gene_set_coverage.csv'))
if(!all(coverage$score_eligible))stop('A set has fewer than 3 expressed genes; inspect coverage before plotting.')
score<-function(E){
 Z<-t(scale(t(E)));Z[!is.finite(Z)]<-0
 zz<-sapply(used,function(g)colMeans(Z[g,,drop=FALSE]))
 ss<-t(gsva(ssgseaParam(E,used,minSize=3,alpha=.25,normalize=TRUE,verbose=FALSE),verbose=FALSE))
 ss<-ss[,names(sets),drop=FALSE]
 list(z_mean=zz,ssgsea=ss)
}
S44<-score(expr44)
# Sensitivity retains the primary gene universe; no additional RA/other samples.
d50<-DGEList(M50[rownames(expr44),,drop=FALSE]);d50<-normLibSizes(d50,method='TMM')
E50<-cpm(d50,log=TRUE,prior.count=1);S50<-score(E50)
fit_models<-function(S,md,cohort){
 ans<-list();k<-1
 for(method in names(S))for(g in names(sets))for(model in c('unadjusted','macrophage','fibroblast','age_sex','mac_age_sex')){
  covars<-switch(model,unadjusted=character(),macrophage='Macrophage',fibroblast='Fibroblast',age_sex=c('age','sex'),mac_age_sex=c('Macrophage','age','sex'))
  # All comparisons involving age use this explicit common complete-case cohort.
  selection<-if(model%in%c('age_sex','mac_age_sex'))complete.cases(md[,c('age','sex')]) else rep(TRUE,nrow(md))
  dat<-md[selection,,drop=FALSE]; yy<-S[[method]][selection,g];sdy<-sd(yy)
  dat$y<-(yy-mean(yy))/sdy
  u<-lm(y~OA,dat);a<-lm(reformulate(c('OA',covars),'y'),dat)
  ci<-confint(a,'OA',level=.95);cf<-summary(a)$coefficients['OA',]
  uc<-summary(u)$coefficients['OA',];uci<-confint(u,'OA')
  X<-model.matrix(a);bread<-solve(crossprod(X));hh<-hatvalues(a)
  meat<-crossprod(X,X*(resid(a)/(1-hh))^2)
  robust_se<-sqrt(diag(bread%*%meat%*%bread))['OA']
  robust_p<-2*pt(-abs(cf[1]/robust_se),df=df.residual(a))
  ans[[k]]<-data.table(cohort=cohort,method=method,geneset=g,model=model,n=nobs(a),n_normal=sum(dat$OA==0),n_OA=sum(dat$OA==1),score_sd=sdy,beta=cf[1],se=cf[2],ci_low=ci[1],ci_high=ci[2],p=cf[4],baseline_beta=uc[1],baseline_ci_low=uci[1],baseline_ci_high=uci[2],baseline_p=uc[4],delta_beta=uc[1]-cf[1],share_if_defined=if(abs(uc[1])>1e-8)(uc[1]-cf[1])/uc[1] else NA_real_)
  ans[[k]][,`:=`(HC3_se=robust_se,HC3_p=robust_p,HC3_ci_low=cf[1]-qt(.975,df.residual(a))*robust_se,HC3_ci_high=cf[1]+qt(.975,df.residual(a))*robust_se,max_leverage=max(hh),max_Cooks_D=max(cooks.distance(a)))]
  k<-k+1
 }
 r<-rbindlist(ans);r[,q:=p.adjust(p,'BH'),by=.(cohort,method,model)];r[,baseline_q:=p.adjust(baseline_p,'BH'),by=.(cohort,method,model)]
 r[,HC3_q:=p.adjust(HC3_p,'BH'),by=.(cohort,method,model)];r
}
fits<-rbind(fit_models(S44,meta44,'primary_QC44'),fit_models(S50,meta50,'sensitivity_all50'))
fwrite(fits,file.path(out,'all_models_standardized.csv'))
# Panel A: all three adjustment factors compared on the SAME complete cases.
cc<-complete.cases(meta44[,c('age','sex')]);mdcc<-meta44[cc,,drop=FALSE]
Scc<-lapply(S44,function(x)x[cc,,drop=FALSE])
Acc<-fit_models(Scc,mdcc,'complete_case_QC')
fwrite(Acc,file.path(out,'complete_case_models.csv'))
A<-Acc[model%in%c('macrophage','fibroblast','age_sex')]
fwrite(A,file.path(out,'Fig2A_all39_coefficient_changes.csv'))
B<-fits[cohort=='primary_QC44' & model%in%c('unadjusted','macrophage') & geneset%in%names4]
fwrite(B,file.path(out,'Fig2B_four_sets_coefficients.csv'))
C<-fits[cohort=='primary_QC44' & model=='macrophage']
C[,display_order:=match(geneset,sort(names(sets)))]
fwrite(C,file.path(out,'Fig2C_all39_adjusted_coefficients.csv'))
summary<-fits[,.(n=n[1],n_normal=n_normal[1],n_OA=n_OA[1],n_tested=.N,n_BH_lt05=sum(q<.05),n_HC3_BH_lt05=sum(HC3_q<.05),n_baseline_BH_lt05=sum(baseline_q<.05),median_delta_all39=median(delta_beta),median_share_baseline_sig=if(any(baseline_q<.05))median(share_if_defined[baseline_q<.05])else NA_real_),by=.(cohort,method,model)]
fwrite(summary,file.path(out,'model_summary.csv'))
Asummary<-A[,.(n=n[1],n_normal=n_normal[1],n_OA=n_OA[1],median_delta=median(delta_beta),Q1=quantile(delta_beta,.25),Q3=quantile(delta_beta,.75)),by=.(method,model)]
fwrite(Asummary,file.path(out,'Fig2A_descriptive_summary.csv'))
meta50$primary_included<-qc
meta50$exclusion_reason<-ifelse(qc,'','deconvolution fit_r < 0.6 or missing')
fwrite(as.data.table(meta50),file.path(out,'sample_inclusion_and_metadata.csv'))
for(nm in names(S44))fwrite(cbind(data.table(sample_id=meta44$sample_id),as.data.table(S44[[nm]])),file.path(out,paste0('scores_primary44_',nm,'.csv')))
fwrite(data.table(sample_id=colnames(expr44),raw_library_size=colSums(M44),filtered_library_size=d44$samples$lib.size,TMM_factor=d44$samples$norm.factors,effective_library_size=d44$samples$lib.size*d44$samples$norm.factors),file.path(out,'normalization_factors_primary44.csv'))
# Persist the analyzed expression matrix and fixed gene-set membership for reuse.
fwrite(cbind(data.table(gene=rownames(expr44)),as.data.table(expr44)),file.path(out,'TMM_logCPM_primary44.csv.gz'))
fwrite(rbindlist(lapply(names(used),function(g)data.table(geneset=g,gene=used[[g]]))),file.path(out,'scored_gene_membership.csv'))
saveRDS(list(scores44=S44,scores50=S50,metadata44=meta44,metadata50=meta50,gene_sets=used,expression44=expr44,models=fits),file.path(out,'Fig2_revision_objects.rds'))
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'))
cat('INPUT',nrow(M),'genes',ncol(M),'samples; duplicated IDs',duplicate_count,'\n')
cat('PRIMARY',table(meta44$disease),'COMPLETE CASES',table(mdcc$disease),'\n')
cat('EXPRESSED_UNIVERSE',nrow(expr44),'SETS',length(used),'COVERAGE_RANGE',range(lengths(used)),'\n')
print(summary);print(Asummary)
stopifnot(nrow(C)==78,nrow(B)==16,nrow(A)==234,all(is.finite(fits$beta)),all(fits$ci_low<=fits$ci_high),all(fits$n==fits$n_normal+fits$n_OA))
cat('VALIDATION_OK\n')
