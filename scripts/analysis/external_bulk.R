Sys.setlocale('LC_ALL','English_United States.utf8')
args<-commandArgs(TRUE); root<-normalizePath(args[1],winslash='/'); out<-file.path(root,'analysis');dir.create(out,showWarnings=FALSE)
suppressPackageStartupMessages({library(data.table);library(edgeR);library(GSVA);library(quadprog)})
M<-Sys.getenv('OA_REPO_ROOT',unset=getwd());N<-file.path(M,'data/external_bulk/source');dir.create(file.path(root,'figures'),recursive=TRUE,showWarnings=FALSE)
readgz<-function(p)fread(text=paste(readLines(gzfile(p),warn=FALSE),collapse='\n'))
genes15<-c('PYCARD','OLR1','CD68','CD14','CSF1R','ITGAM','MRC1','TREM2','AIF1','TYROBP','SPP1','COL1A1','PRG4','THY1','PTPRC')
o<-readRDS(file.path(M,'data/components/mixtures/inputs/donor_count_profiles.rds'))
md<-o$lineage_metadata;lab<-as.character(md$celltype);lab[lab%in%c('T_cell','NK')]<-'T_NK';lab[lab%in%c('B_cell','Plasma','Mast')]<-'B_Plasma_Mast'
ref<-sapply(sort(unique(lab)),function(x)rowSums(o$lineage_sum[,lab==x,drop=FALSE]));rownames(ref)<-rownames(o$lineage_sum)
ref[,'Macrophage']<-rowSums(o$clean_mac_sum);ref<-sweep(ref,2,colSums(ref),'/')*1e6
universe<-readgz(file.path(M,'data/components/bulk/analysis/TMM_logCPM_primary44.csv.gz'))[[1]]
ll<-strsplit(readLines(file.path(M,'data/external_bulk/inputs/macrophage_39.gmt')),'\t');sets<-setNames(lapply(ll,function(x)unique(x[-c(1,2)])),sapply(ll,'[',1))
frozen<-fread(file.path(M,'data/external_bulk/inputs/discovery_scored_gene_membership.csv'));sets<-lapply(names(sets),function(n)frozen[geneset==n,gene])|>setNames(names(sets))
fitcomp<-function(linear,exclude=character()){
 common<-intersect(rownames(ref),rownames(linear));S<-ref[common,,drop=FALSE];B<-linear[common,,drop=FALSE]
 spec<-S/rowMeans(S);spec[!is.finite(spec)]<-0
 markers<-unique(unlist(lapply(colnames(S),function(cl){ix<-which(S[,cl]>10 & !rownames(S)%in%exclude);rownames(S)[ix[order(spec[ix,cl],decreasing=TRUE)[seq_len(min(120,length(ix)))]]]})))
 markers<-markers[is.finite(rowMeans(B[markers,,drop=FALSE])) & rowMeans(B[markers,,drop=FALSE])>0]
 S<-S[markers,,drop=FALSE];B<-B[markers,,drop=FALSE]
 S<-S*(rowMeans(B)/pmax(rowMeans(S),1e-8));S<-sweep(S,2,colSums(S),'/');B<-sweep(B,2,colSums(B),'/')
 mu<-rowMeans(S);sdv<-apply(S,1,sd);sdv[sdv<1e-8]<-1
 SZ<-(S-mu)/sdv;BZ<-(B-mu)/sdv;k<-ncol(S)
 est<-t(apply(BZ,2,function(y)solve.QP(crossprod(SZ)+diag(1e-8,k),as.vector(crossprod(SZ,y)),cbind(rep(1,k),diag(k)),c(1,rep(0,k)),meq=1)$solution));colnames(est)<-colnames(S)
 r<-sapply(seq_len(ncol(B)),function(j)cor(B[,j],S%*%est[j,]))
 list(est=est,fit_r=r,markers=markers)
}
fitone<-function(y,OA,mac,fib,model){
 if(!is.finite(sd(y))||sd(y)<1e-10)return(NULL)
 d<-data.frame(y=as.vector(scale(y)),OA=OA,mac=mac,fib=fib)
 formula<-switch(model,unadjusted=y~OA,macrophage=y~OA+mac,fibroblast=y~OA+fib)
 a<-lm(formula,d);co<-summary(a)$coef
 if(!'OA'%in%rownames(co))return(NULL)
 ci<-confint(a,'OA');X<-model.matrix(a)
 hc<-tryCatch({bread<-solve(crossprod(X));meat<-crossprod(X,X*(resid(a)/(1-hatvalues(a)))^2);sqrt(diag(bread%*%meat%*%bread))['OA']},error=function(e)NA_real_)
 data.table(beta=co['OA',1],se=co['OA',2],ci_low=ci[1],ci_high=ci[2],p=co['OA',4],HC3_se=hc,HC3_p=2*pt(-abs(co['OA',1]/hc),df=df.residual(a)),max_leverage=max(hatvalues(a)),condition_number=kappa(X),n=nobs(a))
}
results<-list();candall<-list();qc<-list();coverall<-list();allcomp<-list();candidatesens<-list()
for(ds in c('GSE82107','GSE143514')){
 if(ds=='GSE82107'){
  tab<-fread(file.path(M,'data/external_bulk/inputs/GSE82107_probe_expression.csv'));probe<-tab[[1]];B<-as.matrix(tab[,-1]);storage.mode(B)<-'double';rownames(B)<-probe
  mp<-as.data.frame(fread(file.path(M,'data/external_bulk/inputs/GPL570_probe_gene_mapping.csv')));mp<-mp[match(probe,mp$probe),,drop=FALSE];good<-!is.na(mp$gene);B<-B[good,,drop=FALSE];gg<-mp$gene[good]
  original_quantiles<-quantile(B,c(0,.01,.5,.99,1));log_applied<-quantile(B,.99)>100
  if(log_applied)B<-log2(pmax(B,1))
  E<-rowsum(B,gg,reorder=FALSE)/as.vector(table(factor(gg,levels=unique(gg))))
  # A phenotype-blind alternative retains the probe with largest overall mean.
  index<-order(rowMeans(B),decreasing=TRUE);pick<-index[!duplicated(gg[index])];Ealt<-B[pick,,drop=FALSE];rownames(Ealt)<-gg[pick]
  OA<-as.integer(as.numeric(sub('GSM','',colnames(E)))>=2183539);linear<-2^E
  norminfo<-list(submitter='dChip invariant-set normalization and model-based expression',log2_applied=log_applied,input_quantiles=original_quantiles)
 }else{
  tab<-readgz(file.path(N,ds,'GSE143514_mRNA_raw_count.txt.gz'));B<-as.matrix(tab[,-1]);storage.mode(B)<-'double';rownames(B)<-tab[[1]]
  stopifnot(all(is.finite(B)),all(B>=0));fractional_entries<-sum(B!=round(B))
  if(anyDuplicated(rownames(B)))B<-rowsum(B,rownames(B),reorder=FALSE)
  OA<-as.integer(grepl('Osteoarthritis',colnames(B)));dge<-DGEList(B);keep<-filterByExpr(dge,group=OA);dge<-dge[keep,,keep.lib.sizes=FALSE];dge<-normLibSizes(dge,method='TMM');E<-cpm(dge,log=TRUE,prior.count=1);linear<-cpm(dge,log=FALSE)
  Ealt<-log2(sweep(B,2,colSums(B),'/')*1e6+1)
  norminfo<-list(input_genes=nrow(B),filtered_genes=nrow(E),fractional_entries=fractional_entries,normalization=dge$samples,source_note='Submitter labels raw counts but reports Cuffquant/Cuffdiff and FPKM in methods; fractional values preserved, no rounding or count-likelihood test applied')
 }
 stopifnot(length(OA)==ncol(E),length(unique(OA))==2,all(is.finite(E)))
 comp<-fitcomp(linear);comp15<-fitcomp(linear,genes15)
 metadata<-data.table(dataset=ds,sample=colnames(E),OA=OA,group=ifelse(OA==1,'OA','Control'),fit_r=comp$fit_r,as.data.table(comp$est))
 metadata15<-data.table(dataset=ds,sample=colnames(E),OA=OA,fit_r=comp15$fit_r,as.data.table(comp15$est))
 fwrite(metadata,file.path(out,paste0(ds,'_composition_metadata.csv')));fwrite(metadata15,file.path(out,paste0(ds,'_composition_exclude15.csv')))
 fwrite(data.table(gene=comp$markers,in_panel15=comp$markers%in%genes15),file.path(out,paste0(ds,'_composition_markers.csv')))
 Ee<-E[intersect(universe,rownames(E)),,drop=FALSE];used<-lapply(sets,intersect,y=rownames(Ee));coverage<-data.table(dataset=ds,geneset=names(sets),n_discovery=lengths(sets),n_external=lengths(used),coverage_fraction=lengths(used)/lengths(sets),eligible=lengths(used)>=3);coverall[[ds]]<-coverage;used<-used[lengths(used)>=3]
 Z<-t(scale(t(Ee)));Z[!is.finite(Z)]<-0
 scores<-list(z_mean=t(sapply(used,function(g)colMeans(Z[g,,drop=FALSE]))),ssgsea=gsva(ssgseaParam(Ee,used,minSize=3,alpha=.25,normalize=TRUE,verbose=FALSE),verbose=FALSE))
 for(method in names(scores)){
  S<-scores[[method]];fwrite(data.table(geneset=rownames(S),as.data.table(S)),file.path(out,paste0(ds,'_scores_',method,'.csv')))
  for(model in c('unadjusted','macrophage','fibroblast'))for(g in rownames(S)){
   a<-fitone(S[g,],OA,metadata$Macrophage,metadata$Fibroblast,model)
   if(!is.null(a)){a[,`:=`(dataset=ds,method=method,geneset=g,model=model)];results[[length(results)+1]]<-a}
   if(model=='macrophage'){
    a<-fitone(S[g,],OA,metadata15$Macrophage,metadata15$Fibroblast,model)
    if(!is.null(a)){a[,`:=`(dataset=ds,method=method,geneset=g,model='macrophage_exclude15')];results[[length(results)+1]]<-a}
   }
  }
 }
 for(g in intersect(genes15,rownames(E)))for(model in c('unadjusted','macrophage','fibroblast')){
  a<-fitone(E[g,],OA,metadata$Macrophage,metadata$Fibroblast,model)
  if(!is.null(a)){a[,`:=`(dataset=ds,gene=g,model=model,rho_mac=cor(E[g,],metadata$Macrophage,method='spearman'))];candall[[length(candall)+1]]<-a}
  a<-fitone(E[g,],OA,metadata15$Macrophage,metadata15$Fibroblast,model)
  if(!is.null(a)){a[,`:=`(dataset=ds,gene=g,model=paste0(model,'_exclude15'))];candidatesens[[length(candidatesens)+1]]<-a}
  if(g%in%rownames(Ealt)){
   a<-fitone(Ealt[g,],OA,metadata$Macrophage,metadata$Fibroblast,model)
   if(!is.null(a)){a[,`:=`(dataset=ds,gene=g,model=paste0(model,'_alternative_expression'))];candidatesens[[length(candidatesens)+1]]<-a}
  }
 }
 for(g in intersect(c('OLR1','PYCARD'),rownames(E)))for(j in seq_along(OA)){
  ix<-setdiff(seq_along(OA),j);a<-fitone(E[g,ix],OA[ix],metadata$Macrophage[ix],metadata$Fibroblast[ix],'macrophage')
  if(!is.null(a)){a[,`:=`(dataset=ds,gene=g,omitted_sample=colnames(E)[j])];allcomp[[length(allcomp)+1]]<-a}
 }
 qc[[ds]]<-data.table(dataset=ds,n=ncol(E),n_OA=sum(OA),n_control=sum(OA==0),genes=nrow(E),fixed_universe_genes=nrow(Ee),sets_tested=length(used),min_fit_r=min(comp$fit_r),max_fit_r=max(comp$fit_r),min_mac=min(metadata$Macrophage),max_mac=max(metadata$Macrophage))
 fwrite(data.table(gene=rownames(E),as.data.table(E)),file.path(out,paste0(ds,'_normalized_expression.csv.gz')))
 saveRDS(list(metadata=metadata,expression=E,composition_exclude15=metadata15,normalization=norminfo),file.path(out,paste0(ds,'_objects.rds')))
 variable<-order(apply(Ee,1,var),decreasing=TRUE)[seq_len(min(2000,nrow(Ee)))];pc<-prcomp(t(Ee[variable,,drop=FALSE]),center=TRUE,scale.=FALSE)
 fwrite(data.table(sample=colnames(E),group=metadata$group,PC1=pc$x[,1],PC2=pc$x[,2]),file.path(out,paste0(ds,'_PCA.csv')))
 fwrite(data.table(sample=colnames(E),median_expression=apply(E,2,median),nonfinite_genes=colSums(!is.finite(E))),file.path(out,paste0(ds,'_sample_qc.csv')))
 cat('COMPLETE',ds,ncol(E),'samples',length(used),'sets\n')
}
r<-rbindlist(results);r[,q:=p.adjust(p,'BH'),by=.(dataset,method,model)];r[,HC3_q:=p.adjust(HC3_p,'BH'),by=.(dataset,method,model)]
disc<-fread(file.path(M,'data/external_bulk/inputs/discovery_models.csv'))[cohort=='primary_QC44'&model%in%c('unadjusted','macrophage','fibroblast'),.(method,geneset,model,discovery_beta=beta,discovery_q=q)]
disc_alt<-copy(disc[model=='macrophage']);disc_alt[,model:='macrophage_exclude15'];disc<-rbind(disc,disc_alt)
r<-merge(r,disc,by=c('method','geneset','model'),all.x=TRUE);r[,same_direction_as_discovery:=sign(beta)==sign(discovery_beta)]
fwrite(r,file.path(out,'gene_set_models.csv'));fwrite(rbindlist(coverall),file.path(out,'gene_set_coverage.csv'));fwrite(rbindlist(qc),file.path(out,'cohort_qc.csv'))
candidates<-rbindlist(candall);candidates[,q:=p.adjust(p,'BH'),by=.(dataset,model)];fwrite(candidates,file.path(out,'candidate_gene_models.csv'));fwrite(rbindlist(candidatesens),file.path(out,'candidate_sensitivity.csv'));fwrite(rbindlist(allcomp),file.path(out,'candidate_leave_one_out.csv'))
summary<-r[,.(n=n[1],sets=.N,BH_lt05=sum(q<.05,na.rm=TRUE),HC3_BH_lt05=sum(HC3_q<.05,na.rm=TRUE),direction_agreement=sum(same_direction_as_discovery,na.rm=TRUE)),by=.(dataset,method,model)];fwrite(summary,file.path(out,'model_summary.csv'))
pdf(file.path(root,'figures/External_validation.pdf'),width=9,height=6,onefile=TRUE)
for(ds in c('GSE82107','GSE143514')){
 a<-candidates[dataset==ds & gene%in%c('OLR1','PYCARD','CD68') & model%in%c('unadjusted','macrophage')];a[,label:=paste(gene,model)];a<-a[order(gene,model)]
 par(mar=c(5,11,3,2));plot(a$beta,seq_len(nrow(a)),xlim=range(c(a$ci_low,a$ci_high)),yaxt='n',ylab='',xlab='Standardized OA coefficient (95% CI)',pch=19,main=ds);axis(2,at=seq_len(nrow(a)),labels=a$label,las=1);segments(a$ci_low,seq_len(nrow(a)),a$ci_high,seq_len(nrow(a)));abline(v=0,lty=2,col='gray50')
 p<-fread(file.path(out,paste0(ds,'_PCA.csv')));par(mar=c(5,5,3,2));plot(p$PC1,p$PC2,pch=19,col=ifelse(p$group=='OA','#D55E00','#0072B2'),xlim=extendrange(p$PC1,f=.08),ylim=extendrange(p$PC2,f=.15),xlab='PC1',ylab='PC2',main=paste(ds,'expression PCA'));text(p$PC1,p$PC2,labels=seq_len(nrow(p)),pos=3,cex=.7);legend('topleft',c('Control','OA'),col=c('#0072B2','#D55E00'),pch=19,bty='n')
}
dev.off();writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'));print(summary);print(candidates[gene%in%c('OLR1','PYCARD'),.(dataset,gene,model,beta,ci_low,ci_high,p,q,rho_mac)])
