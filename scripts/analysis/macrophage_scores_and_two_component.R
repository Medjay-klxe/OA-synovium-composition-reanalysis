args<-commandArgs(trailingOnly=TRUE);inp<-if(length(args)>0)args[1]else'work/fig6_revision/inputs';out<-if(length(args)>1)args[2]else'work/fig6_revision/analysis';dir.create(out,recursive=TRUE,showWarnings=FALSE)
suppressPackageStartupMessages(library(data.table))
gs<-fread(file.path(inp,'gene_set_mitoxyperilysis.csv'));setnames(gs,1,'module')
exec<-c('M2_mito_damage','M3_antiox_buffer','M4_mtorc2','M4_mtorc1','M5_cytoskeleton','M6_membrane_lysis')
sigs<-list(Mac_OLR1=c('OLR1','TREM2','GPNMB','ACKR3','GPRC5A','IL1A'),Resident=c('LYVE1','FOLR2','MRC1','F13A1','IGF1','ABCA6'),Mono_FCN1=c('FCN1','VCAN','S100A12','SERPINB2','PLAC8'),LAM_SPP1=c('SPP1','MMP9','MMP7','CHIT1','FABP4'),AntigenPres=c('CD1C','FCER1A','CLEC10A','IL1R2'),Prolif=c('MKI67','TOP2A','HJURP'))
objects<-lapply(c('GSE216651','GSE248453'),function(n)readRDS(file.path(inp,paste0(n,'_target_counts.rds'))));names(objects)<-c('GSE216651','GSE248453')
balanced_reference<-function(X,md){cw<-1/(uniqueN(md$sample)*as.numeric(table(md$sample)[md$sample]));mu<-drop(X%*%cw);sd<-sqrt(drop((X-mu)^2%*%cw));sd[sd<=1e-12]<-1;list(mean=mu,sd=sd)}
make_scores<-function(X,ref){Z<-sweep(sweep(X,1,ref$mean,'-'),1,ref$sd,'/');mod<-sapply(unique(gs$module),function(m){q<-gs[module==m&gene%in%rownames(Z)];if(nrow(q)==0)return(rep(NA_real_,ncol(Z)));colMeans(Z[q$gene,,drop=FALSE]*q$direction)});S<-cbind(Composite=rowMeans(mod[,exec,drop=FALSE]),Composite_noM6=rowMeans(mod[,setdiff(exec,'M6_membrane_lysis'),drop=FALSE]),mod,sapply(sigs,function(gg)colMeans(Z[intersect(gg,rownames(Z)),,drop=FALSE])));list(S=S,mod=mod)}
prepared<-list()
for(cohort in names(objects)){
 o<-objects[[cohort]];md<-copy(o$metadata);X<-log1p(sweep(o$counts,2,o$library_total,'/')*1e4);ref<-balanced_reference(X,md);z<-make_scores(X,ref);prepared[[cohort]]<-list(md=md,X=X,ref=ref,S=z$S)
 fwrite(data.table(gene=rownames(X),mean=ref$mean,sd=ref$sd),file.path(out,paste0(cohort,'_frozen_donor_balanced_reference.csv')))
 fwrite(cbind(md,as.data.table(z$S)),file.path(out,paste0(cohort,'_cell_scores.csv.gz')))
 dn<-cbind(md[,intersect(c('sample','group','pain'),names(md)),with=FALSE],as.data.table(z$S));dn[,OLR1_logNorm:=X['OLR1',]];dn[,OLR1_positive:=as.numeric(o$counts['OLR1',]>0)];dn[,n_cells:=1]
 donor<-dn[,lapply(.SD,mean),by=intersect(c('sample','group','pain'),names(md)),.SDcols=setdiff(names(dn),c('sample','group','pain','n_cells'))];donor[,n_cells:=as.numeric(table(md$sample)[sample])]
 fwrite(donor,file.path(out,paste0(cohort,'_donor_scores.csv')))
 # Reconstruct archived within-cell scores independently of stored composite columns.
 legacyZ<-t(scale(t(o$data)));legacyZ[!is.finite(legacyZ)]<-0
 legacyMod<-sapply(unique(gs$module),function(m)colMeans(legacyZ[intersect(gs[module==m,gene],rownames(legacyZ)),,drop=FALSE]))
 old<-legacyMod[,'M2_mito_damage']+legacyMod[,'M4_mtorc2']-legacyMod[,'M4_mtorc1']+legacyMod[,'M6_membrane_lysis']-legacyMod[,'M3_antiox_buffer']-legacyMod[,'M5_cytoskeleton']
 if(cohort=='GSE216651'){
  stopifnot(max(abs(old-md$ECI_mac))<1e-9)
  fwrite(data.table(cells=nrow(md),max_legacy_score_difference=max(abs(old-md$ECI_mac))),file.path(out,'audit_legacy_cell_score_identity.csv'))
  prepared[[cohort]]$legacy<-old
 }
 if(cohort=='GSE248453'){
  arch<-fread(file.path(inp,'archived_GSE248453_mac_scores.csv'));ii<-match(md$cell,arch$cell);stopifnot(!anyNA(ii))
  fwrite(data.table(cells=nrow(md),max_legacy_score_difference=max(abs(old-arch$ECI[ii]))),file.path(out,'audit_legacy_pain_score_identity.csv'));stopifnot(max(abs(old-arch$ECI[ii]))<1e-9)
 }
}

# Primary estimand: group mean of donor mean cell scores. Equal donor weight and
# cell fractions within donor define a joint donor/subtype distribution.
p<-prepared$GSE216651;md<-p$md;H<-sort(unique(md[group=='Healthy',sample]));O<-sort(unique(md[group=='OA',sample]));donors<-c(H,O);subs<-sort(unique(md$subtype))
subdon<-cbind(md[,.(sample,group,subtype)],as.data.table(p$S))[,c(list(n_cells=.N),lapply(.SD,mean)),by=.(sample,group,subtype),.SDcols=colnames(p$S)]
subdon[,fraction:=n_cells/sum(n_cells),by=sample];stopifnot(nrow(subdon)==48,min(subdon$n_cells)==6)
fwrite(subdon,file.path(out,'subtype_donor_scores.csv'))
decomp<-function(d,value='Composite',hh=H,oo=O){
 part<-function(ids){rbindlist(lapply(subs,function(k){q<-d[subtype==k];ix<-match(ids,q$sample);w<-q$fraction[ix];e<-q[[value]][ix];data.table(subtype=k,w=mean(w),E=sum(w*e)/sum(w))}))}
 a<-part(hh);b<-part(oo);ans<-data.table(subtype=subs,w_H=a$w,w_OA=b$w,E_H=a$E,E_OA=b$E)
 ans[,composition:=(w_OA-w_H)*(E_H+E_OA)/2];ans[,within_profile:=(w_H+w_OA)*(E_OA-E_H)/2];ans[,total:=composition+within_profile]
 dd<-d[,.(score=sum(fraction*get(value))),by=sample];observed<-mean(dd$score[match(oo,dd$sample)])-mean(dd$score[match(hh,dd$sample)])
 stopifnot(abs(sum(ans$total)-observed)<1e-10);ans
}
dc<-decomp(subdon);fwrite(dc,file.path(out,'Fig6B_decomposition_by_subtype.csv'))
totfun<-function(dc)data.table(composition=sum(dc$composition),within_profile=sum(dc$within_profile),total=sum(dc$total))
tt<-totfun(dc);fwrite(tt,file.path(out,'decomposition_total.csv'))
exact<-function(x,g,positive){ii<-which(g==positive);obs<-mean(x[ii])-mean(x[-ii]);cmb<-combn(seq_along(x),length(ii));null<-apply(cmb,2,function(z)mean(x[z])-mean(x[-z]));c(delta=obs,p=mean(abs(null)>=abs(obs)-1e-12),n_allocations=ncol(cmb))}
substats<-rbindlist(lapply(subs,function(k){d<-subdon[subtype==k];q<-exact(d$Composite,d$group,'OA');data.table(subtype=k,mean_H=mean(d[group=='Healthy',Composite]),mean_OA=mean(d[group=='OA',Composite]),delta=q[1],exact_p=q[2],n_allocations=q[3],all_OA_above_all_H=min(d[group=='OA',Composite])>max(d[group=='Healthy',Composite]))}));substats[,p_BH:=p.adjust(exact_p,'BH')]
fwrite(substats,file.path(out,'Fig6A_subtype_contrasts.csv'))
lod<-rbindlist(lapply(donors,function(d){q<-decomp(subdon,hh=setdiff(H,d),oo=setdiff(O,d));q[,omitted_donor:=d];q}));fwrite(lod,file.path(out,'decomposition_LODO.csv'));fwrite(lod[,.(composition=sum(composition),within_profile=sum(within_profile),total=sum(total)),by=omitted_donor],file.path(out,'decomposition_LODO_totals.csv'))
slod<-rbindlist(lapply(donors,function(d)subdon[sample!=d,.(delta=mean(Composite[group=='OA'])-mean(Composite[group=='Healthy'])),by=subtype][,omitted_donor:=d]));fwrite(slod,file.path(out,'subtype_contrasts_LODO.csv'))
draws<-as.matrix(expand.grid(i=1:3,j=1:3,k=1:3));boot<-rbindlist(lapply(seq_len(27),function(a)rbindlist(lapply(seq_len(27),function(b){v<-totfun(decomp(subdon,hh=H[draws[a,]],oo=O[draws[b,]]));v[,resample_id:=(a-1)*27+b];v}))))
fwrite(boot,file.path(out,'decomposition_donor_bootstrap729.csv'));fwrite(melt(boot,id.vars='resample_id')[,.(Q025=quantile(value,.025),Q975=quantile(value,.975)),by=variable],file.path(out,'decomposition_bootstrap_ranges.csv'))
# Legacy pooled reconstruction, and deliberately distinct donor-balanced target.
legacy<-copy(md);legacy[,score:=p$legacy];legacy[,fraction:=1];ld<-legacy[,.(n_cells=.N,Legacy=mean(score)),by=.(sample,group,subtype)];ld[,fraction:=n_cells/sum(n_cells),by=sample]
fwrite(totfun(decomp(ld,value='Legacy')),file.path(out,'legacy_formula_donor_balanced_decomposition.csv'))
pooled<-legacy[,.(n=.N,E=mean(score)),by=.(group,subtype)];pooled[,w:=n/sum(n),by=group];ph<-pooled[group=='Healthy'][match(subs,subtype)];po<-pooled[group=='OA'][match(subs,subtype)];legacydec<-data.table(subtype=subs,composition=(po$w-ph$w)*(po$E+ph$E)/2,within_profile=(po$w+ph$w)*(po$E-ph$E)/2);legacydec[,total:=composition+within_profile];fwrite(legacydec,file.path(out,'audit_legacy_pooled_decomposition.csv'))
# Sensitivity: no singleton M6, common-support 5 subtypes, equal subtype weights,
# reference re-estimation under LODO.
sens<-list(primary=tt,no_M6=totfun(decomp(subdon,'Composite_noM6')))
major<-subdon[,.(min_cells=min(n_cells)),by=subtype][min_cells>=20,subtype]
major_d<-subdon[subtype%in%major];major_d[,fraction:=fraction/sum(fraction),by=sample]
saved_subs<-subs;subs<-major;sens$major5<-totfun(decomp(major_d));subs<-saved_subs
eq<-copy(subdon);eq[,fraction:=1/length(subs)];sens$equal_subtype_weights<-totfun(decomp(eq))
fwrite(rbindlist(sens,idcol='scenario'),file.path(out,'decomposition_sensitivity.csv'))
ref_lodo<-rbindlist(lapply(donors,function(d){keep<-md$sample!=d;rr<-balanced_reference(p$X[,keep,drop=FALSE],md[keep]);ss<-make_scores(p$X[,keep,drop=FALSE],rr)$S[,'Composite'];sd<-copy(md[keep,.(sample,group,subtype)]);sd[,score:=ss];sd<-sd[,.(n_cells=.N,Composite=mean(score)),by=.(sample,group,subtype)];sd[,fraction:=n_cells/sum(n_cells),by=sample];totfun(decomp(sd,hh=setdiff(H,d),oo=setdiff(O,d)))[,omitted_donor:=d]}))
fwrite(ref_lodo,file.path(out,'decomposition_LODO_reestimated_reference.csv'))

# OLR1-associated expression: within-donor positive-minus-negative comparisons,
# with library-depth quartile matching as a technical-depth sensitivity.
coexpr<-c('TREM2','GPNMB','APOE','MARCO','LYVE1','FOLR2','MRC1','F13A1')
cx<-list();detect<-list();coverage<-list()
for(cohort in names(prepared)){
 pp<-prepared[[cohort]];mm<-pp$md;X<-pp$X;pos<-objects[[cohort]]$counts['OLR1',]>0
 for(d in unique(mm$sample)){
  ix<-which(mm$sample==d);ip<-ix[pos[ix]];ineg<-ix[!pos[ix]];stopifnot(length(ip)>=3,length(ineg)>=3)
  detect[[paste(cohort,d)]]<-data.table(cohort=cohort,sample=d,n_cells=length(ix),OLR1_positive_fraction=mean(pos[ix]),median_library=median(mm$library_total[ix]))
  raw<-rowMeans(X[coexpr,ip,drop=FALSE])-rowMeans(X[coexpr,ineg,drop=FALSE])
  bins<-pmin(4,ceiling(rank(mm$library_total[ix],ties.method='first')/length(ix)*4));valid<-sapply(1:4,function(b){v<-ix[bins==b];sum(pos[v])>=3&&sum(!pos[v])>=3});bs<-(1:4)[valid];nv<-sapply(bs,function(b)sum(bins==b));stopifnot(length(bs)>0)
  adj<-Reduce('+',lapply(seq_along(bs),function(k){v<-ix[bins==bs[k]];(rowMeans(X[coexpr,v[pos[v]],drop=FALSE])-rowMeans(X[coexpr,v[!pos[v]],drop=FALSE]))*nv[k]/sum(nv)}))
  coverage[[paste(cohort,d)]]<-data.table(cohort=cohort,sample=d,valid_quartiles=length(bs),covered_cells=sum(nv),total_cells=length(ix))
  cx[[paste(cohort,d)]]<-rbind(data.table(cohort=cohort,sample=d,gene=coexpr,method='within_donor',delta=raw),data.table(cohort=cohort,sample=d,gene=coexpr,method='depth_quartile_matched',delta=adj))
 }
}
cx<-rbindlist(cx);fwrite(cx,file.path(out,'coexpression_donor_contrasts.csv'));fwrite(rbindlist(detect),file.path(out,'OLR1_detection_by_donor.csv'));fwrite(rbindlist(coverage),file.path(out,'coexpression_depth_matching_coverage.csv'))
cxs<-cx[,.(mean_delta=mean(delta),min_donor=min(delta),max_donor=max(delta),n_donors=.N,positive_donors=sum(delta>0),negative_donors=sum(delta<0),LODO_min=min((sum(delta)-delta)/(.N-1)),LODO_max=max((sum(delta)-delta)/(.N-1))),by=.(cohort,gene,method)];fwrite(cxs,file.path(out,'Fig6C_coexpression_summary.csv'))
agreements<-merge(cxs[cohort=='GSE216651'],cxs[cohort=='GSE248453'],by=c('gene','method'),suffixes=c('_216651','_248453'));agreements[,same_sign:=sign(mean_delta_216651)==sign(mean_delta_248453)];fwrite(agreements,file.path(out,'coexpression_cross_cohort_agreement.csv'))

# Pain is a separate clinical contrast: all eight donors have OA.
pain<-fread(file.path(out,'GSE248453_donor_scores.csv'));variables<-c('Composite',names(sigs),'M1_immune_sensing',exec,'S1_metabolic_stress','OLR1_logNorm','OLR1_positive')
variables<-unique(variables);pain_res<-function(dd){rbindlist(lapply(variables,function(v){q<-exact(dd[[v]],dd$pain,'high');data.table(variable=v,delta=q[1],exact_p=q[2],n_allocations=q[3],mean_high=mean(dd[pain=='high'][[v]]),mean_low=mean(dd[pain=='low'][[v]]))}))}
pr<-pain_res(pain);pr[,p_BH:=p.adjust(exact_p,'BH')];fwrite(pr,file.path(out,'pain_exact_permutation_results.csv'))
pl<-rbindlist(lapply(pain$sample,function(d)pain_res(pain[sample!=d])[,omitted_donor:=d]));fwrite(pl,file.path(out,'pain_LODO.csv'))
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'))
saveRDS(list(subdon=subdon,decomposition=dc,total=tt,substats=substats,coexpression=cx,coexpression_summary=cxs,pain_results=pr),file.path(out,'Fig6_revision_objects.rds'))
cat('PRIMARY_DECOMPOSITION\n');print(tt);print(substats);cat('LODO\n');print(lod[,.(composition=sum(composition),within_profile=sum(within_profile),total=sum(total)),by=omitted_donor]);cat('COEXPRESSION\n');print(cxs);cat('PAIN\n');print(pr);cat('ANALYSIS_OK\n')
