# Independent Fig3 audit and revision. No writes to the original project.
# Rscript code/fig3_revision_analysis.R inputs analysis
args<-commandArgs(trailingOnly=TRUE);inp<-if(length(args)>0)args[1]else'inputs';out<-if(length(args)>1)args[2]else'analysis'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
suppressPackageStartupMessages({library(data.table);library(edgeR)})
txt<-readLines(gzfile(file.path(inp,'GSE89408_counts.txt.gz')),warn=FALSE)
ct<-fread(text=paste(txt,collapse='\n'));rm(txt)
M<-as.matrix(ct[,-1]);rownames(M)<-ct[[1]];storage.mode(M)<-'double';rm(ct)
if(anyDuplicated(rownames(M)))stop('Duplicate identifiers require explicit input review.')
md<-as.data.frame(fread(file.path(inp,'deconv_proportions_218.csv')))
names(md)[names(md)=='sample']<-'sample_id';md<-md[match(colnames(M),md$sample_id),c('sample_id','disease','age','sex','fit_r','Macrophage','Fibroblast')]
stopifnot(!anyNA(md$sample_id),!anyDuplicated(md$sample_id),all(is.finite(M)),all(M>=0))
md$OA<-as.integer(md$disease=='Osteoarthritis')
eligible<-md$disease%in%c('Normal','Osteoarthritis');good<-is.finite(md$fit_r)&md$fit_r>=.6
ix50<-which(eligible);ix44<-which(eligible&good);m50<-md[ix50,];m44<-md[ix44,]
genes<-c('PYCARD','OLR1','CD68','CD14','CSF1R','ITGAM','MRC1','TREM2','AIF1','TYROBP','SPP1','COL1A1','PRG4','THY1','PTPRC')
category<-setNames(c(rep('Candidate',2),rep('Macrophage marker',9),rep('Fibroblast marker',3),'Pan-leukocyte'),genes)
auc<-function(s,y){stopifnot(length(s)==length(y),all(is.finite(s)));n1<-sum(y==1);n0<-sum(y==0);(sum(rank(s,ties.method='average')[y==1])-n1*(n1+1)/2)/(n1*n0)}
auc_ci_old<-function(s,y){set.seed(7);v<-replicate(2000,{i<-sample(seq_along(y),replace=TRUE);if(length(unique(y[i]))<2)NA_real_ else auc(s[i],y[i])});unname(quantile(v,c(.025,.975),na.rm=TRUE))}

# Audit the archived point estimates and their fixed-score bootstrap schedule.
old<-fread(file.path(inp,'archived_P14_candidate_genes.csv'))
raw50<-log2(M[genes,ix50,drop=FALSE]+1)
oldcheck<-rbindlist(lapply(genes,function(g){
 y<-raw50[g,];D<-m50$OA;mac<-m50$Macrophage;fib<-m50$Fibroblast
 u<-lm(y~D);a<-lm(y~D+mac);f<-lm(y~D+fib);rm<-resid(lm(y~mac));rf<-resid(lm(y~fib));ci<-confint(a,'D')
 b0<-coef(u)[2];b1<-coef(a)[2];l0<-glm(D~mac,family=binomial);l1<-glm(D~mac+y,family=binomial)
 cr<-auc_ci_old(y,D);ca<-auc_ci_old(rm,D)
 data.table(gene=g,beta_raw=b0,p_raw=summary(u)$coefficients[2,4],beta_adj_mac=b1,p_adj_mac=summary(a)$coefficients[2,4],beta_adj_lo=ci[1],beta_adj_hi=ci[2],attr_mac=(b0-b1)/b0,attr_fib=(b0-coef(f)[2])/b0,rho_mac=cor(y,mac,method='spearman'),auc_raw=auc(y,D),auc_raw_lo=cr[1],auc_raw_hi=cr[2],auc_resid_mac=auc(rm,D),auc_res_lo=ca[1],auc_res_hi=ca[2],auc_resid_fib=auc(rf,D),p_added_to_mac=anova(l0,l1,test='Chisq')$`Pr(>Chi)`[2],auc_mac_plus_gene=auc(fitted(l1),D),delta_auc=auc(fitted(l1),D)-auc(mac,D))
}))
fwrite(oldcheck,file.path(out,'audit_archived_P14_recomputed.csv'))
joined<-merge(oldcheck,old,by='gene',suffixes=c('_new','_old'));auditdiff<-rbindlist(lapply(setdiff(names(oldcheck),'gene'),function(k)data.table(metric=k,max_abs_difference=max(abs(joined[[paste0(k,'_new')]]-joined[[paste0(k,'_old')]])))))
fwrite(auditdiff,file.path(out,'audit_archived_P14_differences.csv'))
defs<-list('OA vs Normal'=eligible,'OA vs all non-OA'=rep(TRUE,nrow(md)),'OA vs Rheumatoid arthritis'=md$OA==1|grepl('Rheumatoid',md$disease))
oldcontrol<-rbindlist(lapply(names(defs),function(nm){ix<-which(defs[[nm]]);D<-md$OA[ix];rbindlist(lapply(c('PYCARD','OLR1','CD68','Macrophage fraction'),function(g){s<-if(g=='Macrophage fraction')md$Macrophage[ix]else log2(M[g,ix]+1);ci<-auc_ci_old(s,D);data.table(definition=nm,gene=g,n=length(D),n_oa=sum(D),auc=auc(s,D),lo=ci[1],hi=ci[2])}))}))
archcontrol<-fread(file.path(inp,'archived_P14b_control_definition.csv'));j<-merge(oldcontrol,archcontrol,by=c('definition','gene'),suffixes=c('_new','_old'))
fwrite(oldcontrol,file.path(out,'audit_archived_P14b_recomputed.csv'))
control_diff<-max(abs(as.matrix(j[,.(auc_new-auc_old,lo_new-lo_old,hi_new-hi_old)])))
cat('ARCHIVE_P14_MAX_DIFF',max(auditdiff$max_abs_difference),'ARCHIVE_CONTROL_MAX_DIFF',control_diff,'\n')

# Normalize the same primary cohort as Fig2, then freeze its gene universe.
normalise<-function(counts,D){d<-DGEList(counts);keep<-filterByExpr(d,group=D);d<-d[keep,,keep.lib.sizes=FALSE];d<-normLibSizes(d,method='TMM');list(d=d,E=cpm(d,log=TRUE,prior.count=1),genes=rownames(d))}
norm44<-normalise(M[,ix44],m44$OA);E44<-norm44$E
stopifnot(nrow(m44)==44,all(genes%in%rownames(E44)))
norm50<-DGEList(M[rownames(E44),ix50,drop=FALSE]);norm50<-normLibSizes(norm50,method='TMM');E50<-cpm(norm50,log=TRUE,prior.count=1)
preproc<-rbindlist(lapply(c('primary_QC44','sensitivity_all50'),function(cohort){
 ix<-if(cohort=='primary_QC44')ix44 else ix50;D<-md$OA[ix];lib<-colSums(M[,ix])
 Es<-list(raw_log2_counts=log2(M[genes,ix]+1),log2_CPM_plus1=log2(sweep(M[genes,ix],2,lib,'/')*1e6+1),TMM_logCPM=if(cohort=='primary_QC44')E44[genes,]else E50[genes,])
 rbindlist(lapply(names(Es),function(prep){rbindlist(lapply(genes[1:2],function(g){
  x<-as.numeric(Es[[prep]][g,]);z<-as.numeric(scale(x));u<-lm(z~D);a<-lm(z~D+md$Macrophage[ix]);ci0<-confint(u,'D');ci1<-confint(a,'D')
  data.table(cohort=cohort,preprocessing=prep,gene=g,n=length(D),beta=coef(u)[2],ci_low=ci0[1],ci_high=ci0[2],p=summary(u)$coefficients[2,4],beta_adj=coef(a)[2],adj_ci_low=ci1[1],adj_ci_high=ci1[2],p_adj=summary(a)$coefficients[2,4],auc=auc(x,D),rho_mac=cor(x,md$Macrophage[ix],method='spearman'),rho_raw_library_size=cor(x,lib,method='spearman'))
 }))}))
}))
fwrite(preproc,file.path(out,'candidate_preprocessing_sensitivity.csv'))
fwrite(preproc[cohort=='primary_QC44'],file.path(out,'Fig3A_preprocessing_coefficients.csv'))
fwrite(cbind(data.table(sample_id=m44$sample_id),as.data.table(t(E44[genes,]))),file.path(out,'primary44_normalized_gene_expression.csv'))
gene_models<-function(E,meta,cohort){rbindlist(lapply(genes,function(g){
  x<-as.numeric(E[g,]);z<-(x-mean(x))/sd(x);dd<-data.frame(y=z,OA=meta$OA,mac=meta$Macrophage,fib=meta$Fibroblast)
  rbindlist(lapply(c('unadjusted','macrophage','fibroblast'),function(model){
   fit<-lm(switch(model,unadjusted=y~OA,macrophage=y~OA+mac,fibroblast=y~OA+fib),dd);cf<-summary(fit)$coefficients['OA',];ci<-confint(fit,'OA')
   X<-model.matrix(fit);bread<-solve(crossprod(X));vv<-bread%*%crossprod(X,X*(resid(fit)/(1-hatvalues(fit)))^2)%*%bread;hc<-sqrt(vv['OA','OA']);hp<-2*pt(-abs(cf[1]/hc),df=df.residual(fit))
   data.table(cohort=cohort,gene=g,category=category[g],model=model,n=nobs(fit),score_sd=sd(x),beta=cf[1],se=cf[2],ci_low=ci[1],ci_high=ci[2],p=cf[4],HC3_se=hc,HC3_p=hp,rho_mac=cor(x,meta$Macrophage,method='spearman'))
  }))
 }))}
fits<-rbind(gene_models(E44,m44,'primary_QC44'),gene_models(E50,m50,'sensitivity_all50'))
fits[,q_panel15:=p.adjust(p,'BH'),by=.(cohort,model)];fits[,HC3_q_panel15:=p.adjust(HC3_p,'BH'),by=.(cohort,model)]
fits[,q_candidates2:=NA_real_];fits[gene%in%genes[1:2],q_candidates2:=p.adjust(p,'BH'),by=.(cohort,model)]
fwrite(fits,file.path(out,'gene_models_all15.csv'))
age_sensitivity<-rbindlist(lapply(c('primary_QC44','sensitivity_all50'),function(cohort){
 meta<-if(cohort=='primary_QC44')m44 else m50;E<-if(cohort=='primary_QC44')E44 else E50
 cc<-complete.cases(meta[,c('age','sex')]);dat<-meta[cc,];dat$sex<-factor(dat$sex)
 rbindlist(lapply(genes[1:2],function(g){dat$y<-as.numeric(scale(E[g,cc]));rbindlist(lapply(c('unadjusted','macrophage','age_sex','mac_age_sex'),function(model){
  rhs<-switch(model,unadjusted='OA',macrophage='OA+Macrophage',age_sex='OA+age+sex',mac_age_sex='OA+Macrophage+age+sex');fit<-lm(as.formula(paste('y~',rhs)),dat);cf<-summary(fit)$coefficients['OA',];ci<-confint(fit,'OA')
  data.table(cohort=cohort,gene=g,model=model,n=nobs(fit),n_normal=sum(dat$OA==0),n_OA=sum(dat$OA==1),beta=cf[1],ci_low=ci[1],ci_high=ci[2],p=cf[4])
 }))}))
}));fwrite(age_sensitivity,file.path(out,'candidate_age_sex_sensitivity.csv'))
fwrite(fits[cohort=='primary_QC44' & model%in%c('unadjusted','macrophage')],file.path(out,'benchmark_all15_coefficients.csv'))
gene_summary<-merge(fits[model=='unadjusted',.(cohort,gene,beta_raw=beta,p_raw=p,q_raw=q_panel15)],fits[model=='macrophage',.(cohort,gene,beta_adj=beta,p_adj=p,q_adj=q_panel15,rho_mac)],by=c('cohort','gene'))
gene_summary[,delta_beta:=beta_raw-beta_adj];gene_summary[,descriptive_share:=(beta_raw-beta_adj)/beta_raw]
fwrite(gene_summary,file.path(out,'coefficient_changes.csv'))

# Reference-anchored projection: pairwise TMM ratio cancels the pair's geometric-mean scaling.
# The reference and its factor are learned ONLY in the training/reference cohort.
project<-function(train_d,test_counts,ref){
 test_counts<-test_counts[rownames(train_d),,drop=FALSE];ref_counts<-train_d$counts[,ref]
 test_lib<-colSums(test_counts);test_factor<-vapply(seq_len(ncol(test_counts)),function(j){pair<-DGEList(cbind(ref_counts,test_counts[,j]));pair<-normLibSizes(pair,method='TMM',refColumn=1);train_d$samples$norm.factors[ref]*(pair$samples$norm.factors[2]/pair$samples$norm.factors[1])},numeric(1))
 all_eff<-c(train_d$samples$lib.size*train_d$samples$norm.factors,test_lib*test_factor)
 # Match edgeR's prior-count scaling to the TRAINING mean effective library size.
 prior_scaled<-all_eff/mean(train_d$samples$lib.size*train_d$samples$norm.factors)
 eff_test<-tail(all_eff,ncol(test_counts));prior_test<-tail(prior_scaled,ncol(test_counts))
 ans<-log2(sweep(sweep(test_counts,2,prior_test,'+'),2,eff_test+2*prior_test,'/')*1e6)
 list(E=ans,factor=test_factor,library_size=test_lib,reference=ref)
}
# Re-express primary training scores using the same formula as projection (exact edgeR check).
f75<-getFromNamespace('.calcFactorQuantile','edgeR')(norm44$d$counts,lib.size=norm44$d$samples$lib.size,p=.75)
ref44<-if(mean(f75)<1e-20)which.max(colSums(sqrt(norm44$d$counts)))else which.min(abs(f75-mean(f75)))
normcheck<-normLibSizes(norm44$d,method='TMM',refColumn=ref44)
stopifnot(max(abs(normcheck$samples$norm.factors-norm44$d$samples$norm.factors))<1e-12)
eff<-norm44$d$samples$lib.size*norm44$d$samples$norm.factors;pc<-eff/mean(eff)
manual_logCPM<-log2(sweep(sweep(norm44$d$counts,2,pc,'+'),2,eff+2*pc,'/')*1e6)
stopifnot(max(abs(manual_logCPM-E44))<1e-10)
extra<-setdiff(seq_len(ncol(M)),ix44);pr<-project(norm44$d,M[,extra,drop=FALSE],ref44)
Eall<-matrix(NA_real_,nrow=length(genes),ncol=ncol(M),dimnames=list(genes,colnames(M)))
Eall[,ix44]<-E44[genes,];Eall[,extra]<-pr$E[genes,]
fwrite(cbind(data.table(sample_id=md$sample_id,disease=md$disease,fit_r=md$fit_r,primary_included=seq_len(nrow(md))%in%ix44,QC_included=good),as.data.table(t(Eall))),file.path(out,'normalized_gene_expression_all218.csv'))

# Refit residualization in every stratified bootstrap sample; descriptive AUC only.
set.seed(20261110);nboot<-2000L;strat_ids<-function(D)c(sample(which(D==0),sum(D==0),replace=TRUE),sample(which(D==1),sum(D==1),replace=TRUE))
bootids<-replicate(nboot,strat_ids(m44$OA),simplify=FALSE)
descriptive<-rbindlist(lapply(genes,function(g){y<-as.numeric(E44[g,]);D<-m44$OA;mac<-m44$Macrophage;rr<-resid(lm(y~mac));a0<-auc(y,D);a1<-auc(rr,D);bb<-vapply(bootids,function(i){q<-resid(lm(y[i]~mac[i]));c(raw=auc(y[i],D[i]),residual=auc(q,D[i]))},numeric(2));rbindlist(lapply(1:2,function(k){ci<-quantile(bb[k,],c(.025,.975));data.table(gene=g,transformation=c('normalized_expression','macrophage_residual')[k],auc=c(a0,a1)[k],lo=ci[1],hi=ci[2],bootstrap='stratified; residualization refitted per resample',direction='higher score predicts OA')}))}))
fwrite(descriptive,file.path(out,'descriptive_raw_residual_AUC.csv'))

# Panel C: fixed direction, fixed 44-sample normalization reference, changing controls.
controlnew<-rbindlist(lapply(names(defs),function(nm){ix<-which(defs[[nm]]&good);D<-md$OA[ix];set.seed(20261111);ii<-replicate(nboot,strat_ids(D),simplify=FALSE)
 rbindlist(lapply(c('PYCARD','OLR1','CD68','Macrophage fraction'),function(g){s<-if(g=='Macrophage fraction')md$Macrophage[ix]else Eall[g,ix];v<-vapply(ii,function(i)auc(s[i],D[i]),numeric(1));ci<-quantile(v,c(.025,.975));data.table(definition=nm,gene=g,n=length(ix),n_OA=sum(D),n_control=sum(1-D),auc=auc(s,D),lo=ci[1],hi=ci[2],direction='higher score predicts OA; never auto-flipped')}))
}))
fwrite(controlnew,file.path(out,'Fig3C_control_definition.csv'))

# Fixed-penalty ridge logistic regression. Intercept unpenalized; no tuning on outcomes.
ridgefit<-function(X,D,lambda=1){X<-as.matrix(X);mu<-colMeans(X);sg<-apply(X,2,sd);sg[!is.finite(sg)|sg<1e-8]<-1;Z<-cbind(1,sweep(sweep(X,2,mu,'-'),2,sg,'/'));fn<-function(b)sum(pmax(Z%*%b,0)+log1p(exp(-abs(Z%*%b)))-D*(Z%*%b))+lambda*sum(b[-1]^2)/2;gr<-function(b)drop(crossprod(Z,plogis(Z%*%b)-D))+c(0,lambda*b[-1]);fit<-optim(c(qlogis(mean(D)),rep(0,ncol(X))),fn,gr,method='BFGS',control=list(reltol=1e-10,maxit=500));if(fit$convergence!=0)stop('Ridge optimizer did not converge');list(beta=fit$par,mu=mu,sd=sg)}
ridgepred<-function(fit,X){Z<-cbind(1,sweep(sweep(as.matrix(X),2,fit$mu,'-'),2,fit$sd,'/'));drop(plogis(Z%*%fit$beta))}
model_names<-c('Macrophage fraction','PYCARD alone','OLR1 alone','Macrophage + PYCARD','Macrophage + OLR1')
predictors<-function(E,mac){list(matrix(mac,ncol=1),matrix(E['PYCARD',],ncol=1),matrix(E['OLR1',],ncol=1),cbind(mac,E['PYCARD',]),cbind(mac,E['OLR1',]))}
repeats<-20L;kfold<-5L;lambdas<-c(.1,1,10);D<-m44$OA
pred<-array(NA_real_,dim=c(44,repeats,5,3),dimnames=list(m44$sample_id,paste0('repeat',seq_len(repeats)),model_names,as.character(lambdas)))
foldrows<-list();nfilter<-numeric();cv_factors<-list();ff<-1
cache_path<-file.path(out,'CV_predictions_cache.rds')
cache_key<-list(input_md5=unname(tools::md5sum(file.path(inp,c('GSE89408_counts.txt.gz','deconv_proportions_218.csv')))),repeats=repeats,kfold=kfold,lambdas=lambdas,version='frozen_reference_cv_v1')
resume<-length(args)>=3 && args[3]=='--resume-cv'
if(resume){cache<-readRDS(cache_path);stopifnot(identical(cache$key,cache_key));pred<-cache$pred;foldrows<-cache$foldrows;cv_factors<-cache$cv_factors}else{
for(rep in seq_len(repeats)){
 set.seed(20261112+rep);fold<-integer(length(D));for(c in 0:1){ix<-which(D==c);fold[ix]<-sample(rep(seq_len(kfold),length.out=length(ix)))}
 for(k in seq_len(kfold)){
  te<-which(fold==k);tr<-which(fold!=k)
  nd<-DGEList(M[,ix44[tr],drop=FALSE]);kk<-filterByExpr(nd,group=D[tr]);nd<-nd[kk,,keep.lib.sizes=FALSE]
  ref<-which.min(abs(log(nd$samples$lib.size)-median(log(nd$samples$lib.size))));nd<-normLibSizes(nd,method='TMM',refColumn=ref);Etr<-cpm(nd,log=TRUE,prior.count=1)
  pp<-project(nd,M[,ix44[te],drop=FALSE],ref);Ete<-pp$E;stopifnot(all(genes[1:2]%in%rownames(Etr)))
  Xtr<-predictors(Etr,m44$Macrophage[tr]);Xte<-predictors(Ete,m44$Macrophage[te])
  for(li in seq_along(lambdas))for(mi in seq_along(model_names)){fit<-ridgefit(Xtr[[mi]],D[tr],lambdas[li]);pred[te,rep,mi,li]<-ridgepred(fit,Xte[[mi]])}
  cv_factors[[ff]]<-data.table(repeat_id=rep,fold=k,sample_id=m44$sample_id[te],n_training=length(tr),n_training_genes=nrow(nd),training_reference=m44$sample_id[tr[ref]],test_TMM_factor=pp$factor,test_library_size=pp$library_size)
  ff<-ff+1
 }
 foldrows[[rep]]<-data.table(repeat_id=rep,sample_id=m44$sample_id,OA=D,fold=fold)
 cat('CV_REPEAT',rep,'OF',repeats,'\n');flush.console()
}
saveRDS(list(key=cache_key,pred=pred,foldrows=foldrows,cv_factors=cv_factors),cache_path)
}
stopifnot(all(is.finite(pred)))
fwrite(rbindlist(foldrows),file.path(out,'CV_fold_assignment.csv'));fwrite(rbindlist(cv_factors),file.path(out,'CV_training_and_projection.csv'))
# Export a full training/test fold to permit an independent optimizer/prediction check.
te<-which(foldrows[[1]]$fold==1);tr<-which(foldrows[[1]]$fold!=1)
nd<-DGEList(M[,ix44[tr],drop=FALSE]);kk<-filterByExpr(nd,group=D[tr]);nd<-nd[kk,,keep.lib.sizes=FALSE]
ref<-which.min(abs(log(nd$samples$lib.size)-median(log(nd$samples$lib.size))));nd<-normLibSizes(nd,method='TMM',refColumn=ref)
Et<-cpm(nd,log=TRUE,prior.count=1);Ev<-project(nd,M[,ix44[te],drop=FALSE],ref)$E
xf<-matrix(NA_real_,nrow=44,ncol=2,dimnames=list(m44$sample_id,genes[1:2]));xf[tr,]<-t(Et[genes[1:2],]);xf[te,]<-t(Ev[genes[1:2],])
fwrite(cbind(data.table(sample_id=m44$sample_id,OA=D,split=ifelse(seq_along(D)%in%te,'test','training'),Macrophage=m44$Macrophage),as.data.table(xf)),file.path(out,'CV_fold1_features_for_independent_check.csv'))
cvrep<-rbindlist(lapply(seq_along(lambdas),function(li){
 rbindlist(lapply(seq_along(model_names),function(mi){
  data.table(lambda=lambdas[li],model=model_names[mi],repeat_id=seq_len(repeats),auc=vapply(seq_len(repeats),function(j)auc(pred[,j,mi,li],D),numeric(1)))
 }))
}))
fwrite(cvrep,file.path(out,'CV_AUC_by_repeat.csv'))
primary_pred<-pred[,,,which(lambdas==1)]
set.seed(20261113);ii<-replicate(nboot,strat_ids(D),simplify=FALSE)
bootcv<-vapply(ii,function(i)vapply(seq_along(model_names),function(mi)mean(vapply(seq_len(repeats),function(rep)auc(primary_pred[i,rep,mi],D[i]),numeric(1))),numeric(1)),numeric(length(model_names)))
cvsummary<-rbindlist(lapply(seq_along(model_names),function(mi){v<-cvrep[lambda==1 & model==model_names[mi],auc];ci<-quantile(bootcv[mi,],c(.025,.975));data.table(model=model_names[mi],n=44,lambda=1,mean_CV_AUC=mean(v),repeat_sd=sd(v),repeat_Q1=quantile(v,.25),repeat_Q3=quantile(v,.75),conditional_lo=ci[1],conditional_hi=ci[2],interval_type='stratified participant bootstrap of fixed out-of-fold predictions; training uncertainty not fully captured')}))
fwrite(cvsummary,file.path(out,'Fig3B_CV_AUC.csv'))
delta<-rbindlist(lapply(1:2,function(g){mi<-g+3;v<-cvrep[lambda==1&model==model_names[mi],auc]-cvrep[lambda==1&model==model_names[1],auc];ci<-quantile(bootcv[mi,]-bootcv[1,],c(.025,.975));data.table(gene=genes[g],mean_delta_AUC=mean(v),conditional_lo=ci[1],conditional_hi=ci[2],repeat_sd=sd(v),comparison='Macrophage + gene minus Macrophage alone',interval_type='paired conditional bootstrap of fixed out-of-fold predictions')}))
fwrite(delta,file.path(out,'CV_incremental_AUC_candidates.csv'))
fwrite(cvrep[,.(mean_CV_AUC=mean(auc),repeat_sd=sd(auc)),by=.(lambda,model)],file.path(out,'CV_penalty_sensitivity.csv'))
predlong<-rbindlist(lapply(seq_len(repeats),function(j)rbindlist(lapply(seq_along(model_names),function(mi)data.table(sample_id=m44$sample_id,OA=D,repeat_id=j,model=model_names[mi],prediction=primary_pred[,j,mi])))))
predlong[,prediction_17digits:=sprintf('%.17g',prediction)]
fwrite(predlong,file.path(out,'CV_out_of_fold_predictions.csv'))
fwrite(as.data.table(md[,c('sample_id','disease','age','sex','fit_r','Macrophage','Fibroblast','OA')]),file.path(out,'all_metadata218.csv'))
fwrite(cbind(data.table(sample_id=m44$sample_id),as.data.table(t(E44[genes,]))),file.path(out,'primary44_normalized_gene_expression.csv'))
fwrite(data.table(sample_id=m44$sample_id,filtered_library_size=norm44$d$samples$lib.size,TMM_factor=norm44$d$samples$norm.factors),file.path(out,'primary44_normalization.csv'))
fwrite(data.table(gene=rownames(E44)),file.path(out,'primary_expression_universe.csv'))
saveRDS(list(E44=E44[genes,],E50=E50[genes,],Eall=Eall,md=md,ix44=ix44,ix50=ix50,models=fits,CV_predictions=pred,CV_repeat=cvrep),file.path(out,'Fig3_revision_objects.rds'))
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'))
cat('PRIMARY_GROUPS',table(m44$disease),'ALL_QC_GROUPS\n');print(table(md$disease[good]));cat('CANDIDATE_EFFECTS\n');print(fits[cohort=='primary_QC44'&gene%in%genes[1:2]]);cat('CV_SUMMARY\n');print(cvsummary);print(delta);print(controlnew);cat('DESCRIPTIVE_CANDIDATE_AUC\n');print(descriptive[gene%in%genes[1:2]]);cat('VALIDATION_OK\n')
