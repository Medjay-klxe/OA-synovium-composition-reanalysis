# Fig5 donor-balanced controlled mixtures. Run: Rscript code/fig5_revision_analysis.R inputs analysis
args<-commandArgs(trailingOnly=TRUE);inp<-if(length(args)>0)args[1]else'inputs';out<-if(length(args)>1)args[2]else'analysis';dir.create(out,recursive=TRUE,showWarnings=FALSE)
suppressPackageStartupMessages({library(data.table);library(edgeR)})
o<-readRDS(file.path(inp,'donor_count_profiles.rds'));ref<-fread(file.path(inp,'bulk44_score_reference.csv'));gset<-fread(file.path(inp,'gene_set_mitoxyperilysis.csv'));setnames(gset,1,'module')
don<-as.data.table(o$donor_metadata);H<-don[group=='Healthy',sample];O<-don[group=='OA',sample];donors<-c(H,O);stopifnot(length(H)==3,length(O)==3)
genes<-ref$gene;stopifnot(all(genes%in%rownames(o$lineage_sum)))
anchors<-fread(file.path(inp,'fraction_anchors.csv'));flo<-anchors[scenario=='primary_bulk44',f_low];fhi<-anchors[scenario=='primary_bulk44',f_high]
exec<-c('M2_mito_damage','M3_antiox_buffer','M4_mtorc2','M4_mtorc1','M5_cytoskeleton','M6_membrane_lysis')
candidate<-c('CD68','OLR1','PYCARD','NINJ1','COL1A1')

# Audit exact archived arithmetic; no interpretation is borrowed from legacy labels.
legacy<-readRDS(file.path(inp,'archived_P16b_result.rds'));arch<-fread(file.path(inp,'archived_P16b_factorial.csv'))
audit<-rbindlist(lapply(names(legacy$scores),function(v){x<-legacy$scores[[v]];dd<-legacy$design;f<-sort(unique(dd$f));hl<-mean(x[dd$state=='H'&dd$f==f[1]]);hh<-mean(x[dd$state=='H'&dd$f==f[2]]);ol<-mean(x[dd$state=='OA'&dd$f==f[1]]);oh<-mean(x[dd$state=='OA'&dd$f==f[2]]);data.table(variable=v,mean_fLO_H=hl,mean_fHI_H=hh,mean_fLO_OA=ol,mean_fHI_OA=oh,d_comp=(hh+oh-hl-ol)/2,d_state=(ol+oh-hl-hh)/2,d_inter=oh-ol-hh+hl)}))
aa<-merge(audit,arch,by='variable',suffixes=c('_new','_archive'));dif<-rbindlist(lapply(setdiff(names(audit),'variable'),function(k)data.table(metric=k,max_abs_difference=max(abs(aa[[paste0(k,'_new')]]-aa[[paste0(k,'_archive')]])))))
fwrite(audit,file.path(out,'audit_archived_factorial_recomputed.csv'));fwrite(dif,file.path(out,'audit_archived_factorial_differences.csv'))
fwrite(as.data.table(table(legacy$design$state,legacy$design$f)),file.path(out,'audit_legacy_simulation_counts.csv'))
lc<-as.data.table(legacy$curve);lcm<-melt(lc,id.vars='f',variable.name='variable',value.name='value');fwrite(lcm[,.(mean=mean(value),sd=sd(value),n=.N),by=.(variable,f)],file.path(out,'audit_archived_curve_recomputed.csv'))

# Mean raw UMI count per recovered cell, separately by donor and compartment/subtype.
P<-sweep(o$subtype_sum,2,o$subtype_n[colnames(o$subtype_sum)],'/')
rest<-sapply(donors,function(d){ix<-which(o$lineage_metadata$sample==d&o$lineage_metadata$celltype!='Macrophage');rowSums(o$lineage_sum[,ix,drop=FALSE])/sum(o$lineage_metadata$n_cells[ix])})
colnames(rest)<-donors
mac_natural<-sweep(o$clean_mac_sum[,donors,drop=FALSE],2,o$clean_mac_n[donors],'/')
subtypes<-sort(unique(o$subtype_metadata$subtype));ntab<-matrix(0,nrow=length(subtypes),ncol=6,dimnames=list(subtypes,donors))
for(i in seq_len(nrow(o$subtype_metadata))){rr<-o$subtype_metadata[i];ntab[rr$subtype,rr$sample]<-rr$n_cells}
stopifnot(all(ntab>0),identical(as.numeric(colSums(ntab)),as.numeric(o$clean_mac_n[donors])))
w<-rowMeans(sweep(ntab,2,colSums(ntab),'/'));major<-rownames(ntab)[apply(ntab,1,min)>=20];wm<-w[major]/sum(w[major]);stopifnot(length(major)>=3)
make_standardized<-function(weights,equalRNA=FALSE){sapply(donors,function(d){ids<-paste(d,names(weights),sep='::');pp<-P[,ids,drop=FALSE];if(equalRNA)pp<-sweep(pp,2,colSums(pp),'/')*10000;drop(pp%*%weights)})}
mac_std<-make_standardized(w);mac_major<-make_standardized(wm);mac_equalRNA<-make_standardized(w,TRUE);rest_equalRNA<-sweep(rest,2,colSums(rest),'/')*10000
colnames(mac_std)<-colnames(mac_major)<-colnames(mac_equalRNA)<-donors
fwrite(data.table(subtype=names(w),weight=w,minimum_cells_per_donor=apply(ntab,1,min),major_common_support=names(w)%in%major),file.path(out,'fixed_subtype_weights.csv'))
fwrite(as.data.table(ntab,keep.rownames='subtype'),file.path(out,'clean_mac_subtype_counts.csv'))
fwrite(data.table(donor=donors,group=don$group[match(donors,don$sample)],natural_mac_UMI_per_cell=colSums(mac_natural),standardized_mac_UMI_per_cell=colSums(mac_std),rest_UMI_per_cell=colSums(rest)),file.path(out,'donor_RNA_content_profiles.csv'))
# A target+OTHER representation exactly preserves full-library denominators for CPM.
small<-function(M){stopifnot(all(genes%in%rownames(M)));ans<-rbind(M[genes,,drop=FALSE],OTHER=colSums(M)-colSums(M[genes,,drop=FALSE]));stopifnot(all(ans>=-1e-8));ans}
Ms<-small(mac_std);Rs<-small(rest);Mns<-small(mac_natural);Mmaj<-small(mac_major);Meq<-small(mac_equalRNA);Req<-small(rest_equalRNA)
weights<-gset[module%in%c(exec,'M1_immune_sensing'),.(module,gene,direction)];stopifnot(all(weights$gene%in%genes))
fwrite(weights,file.path(out,'scoring_module_directions.csv'))
score_log<-function(E){
 Z<-sweep(sweep(E[genes,,drop=FALSE],1,ref$bulk44_mean,'-'),1,ref$bulk44_sd,'/')
 ms<-do.call(cbind,lapply(exec,function(m){q<-weights[module==m];colMeans(Z[q$gene,,drop=FALSE]*q$direction)}));colnames(ms)<-exec
 rr<-weights[module=='M1_immune_sensing'];immune<-colMeans(Z[rr$gene,,drop=FALSE]*rr$direction)
 S<-cbind(Composite=rowMeans(ms),t(Z[candidate,,drop=FALSE]),M1_immune_sensing=immune)
 list(scores=S,modules=ms,expression=E[candidate,,drop=FALSE])
}
logcpm<-function(C)log2(sweep(C,2,colSums(C),'/')*1e6+1)
cell_counts<-function(mh,mo,bg,lo=flo,hi=fhi){N<-4000;C<-N*cbind(H_low=(1-lo)*bg+lo*mh,H_high=(1-hi)*bg+hi*mh,OA_low=(1-lo)*bg+lo*mo,OA_high=(1-hi)*bg+hi*mo);C}
effects<-function(S){data.table(variable=colnames(S),composition=(S[2,]+S[4,]-S[1,]-S[3,])/2,profile=(S[3,]+S[4,]-S[1,]-S[2,])/2,interaction=S[4,]-S[3,]-S[2,]+S[1,],total=S[4,]-S[1,],composition_at_H=S[2,]-S[1,],composition_at_OA=S[4,]-S[3,],profile_at_low=S[3,]-S[1,],profile_at_high=S[4,]-S[2,])}
evaluate<-function(Mac,Rest,hh=H,oo=O,lo=flo,hi=fhi){C<-cell_counts(rowMeans(Mac[,hh,drop=FALSE]),rowMeans(Mac[,oo,drop=FALSE]),rowMeans(Rest[,hh,drop=FALSE]),lo,hi);sc<-score_log(logcpm(C));list(counts=C,score=sc,effects=effects(sc$scores))}
primary<-evaluate(Ms,Rs);ef<-primary$effects
fwrite(ef,file.path(out,'factorial_primary_effects.csv'))
corners<-as.data.table(primary$score$scores,keep.rownames='condition');fwrite(corners,file.path(out,'factorial_primary_corners.csv'))
stopifnot(max(abs(ef$composition+ef$profile-ef$total))<1e-12)

# Exact empirical donor bootstrap: 3^3 ordered draws in each group, paired healthy background.
draws<-as.matrix(expand.grid(i=1:3,j=1:3,k=1:3));boot<-vector('list',nrow(draws)^2);kk<-1
for(a in seq_len(nrow(draws)))for(b in seq_len(nrow(draws))){z<-evaluate(Ms,Rs,H[draws[a,]],O[draws[b,]])$effects;z[,bootstrap_id:=kk];boot[[kk]]<-z;kk<-kk+1}
boot<-rbindlist(boot);fwrite(boot,file.path(out,'donor_bootstrap_effects_729.csv'))
long<-melt(ef,id.vars='variable',measure.vars=c('composition','profile','interaction','total'),variable.name='effect_type',value.name='effect')
bl<-melt(boot,id.vars=c('variable','bootstrap_id'),measure.vars=c('composition','profile','interaction','total'),variable.name='effect_type',value.name='effect')
bi<-bl[,.(bootstrap_Q025=quantile(effect,.025),bootstrap_Q975=quantile(effect,.975)),by=.(variable,effect_type)]
long<-merge(long,bi,by=c('variable','effect_type'))
lodo<-rbindlist(lapply(donors,function(d){hh<-setdiff(H,d);oo<-setdiff(O,d);z<-evaluate(Ms,Rs,hh,oo)$effects;z[,omitted_donor:=d];z}))
fwrite(lodo,file.path(out,'leave_one_donor_out_effects.csv'))
ll<-melt(lodo,id.vars=c('variable','omitted_donor'),measure.vars=c('composition','profile','interaction','total'),variable.name='effect_type',value.name='effect');lr<-ll[,.(LODO_min=min(effect),LODO_max=max(effect),LODO_positive=sum(effect>0),LODO_negative=sum(effect<0)),by=.(variable,effect_type)]
long<-merge(long,lr,by=c('variable','effect_type'));fwrite(long,file.path(out,'effect_summary_with_sensitivity.csv'));fwrite(long[variable%in%c('Composite','CD68','OLR1','PYCARD')&effect_type%in%c('composition','profile')],file.path(out,'Fig5A_factorial_effects.csv'))

# Donor-level direct macrophage comparison with fixed eight-subtype weights.
don_sc<-score_log(logcpm(Ms));donexpr<-rbindlist(lapply(seq_along(donors),function(j){data.table(donor=donors[j],group=ifelse(donors[j]%in%H,'Healthy','OA'),variable=colnames(don_sc$scores),standardized_score=don_sc$scores[j,],logCPM=c(NA_real_,as.numeric(don_sc$expression[,j]),NA_real_))}))
fwrite(donexpr,file.path(out,'fixed_subtype_donor_scores.csv'));fwrite(donexpr[variable=='OLR1'],file.path(out,'Fig5B_donor_OLR1.csv'))
perms<-combn(seq_len(6),3);pex<-rbindlist(lapply(colnames(don_sc$scores),function(v){y<-don_sc$scores[,v];observed<-mean(y[4:6])-mean(y[1:3]);pv<-apply(perms,2,function(i)mean(y[i])-mean(y[-i]));data.table(variable=v,donor_group_difference=observed,exact_two_sided_permutation_p=mean(abs(pv)>=abs(observed)-1e-12),n_label_allocations=ncol(perms),all_OA_above_all_H=min(y[4:6])>max(y[1:3]))}))
fwrite(pex,file.path(out,'donor_group_exact_permutation_checks.csv'))
# Report raw within-subtype donor expression, not a pooled-cell significance test.
subexpr<-logcpm(o$subtype_sum);rows<-rbindlist(lapply(seq_len(nrow(o$subtype_metadata)),function(j){m<-o$subtype_metadata[j];data.table(donor=m$sample,group=m$group,subtype=m$subtype,n_cells=m$n_cells,gene=candidate,logCPM=subexpr[candidate,j])}))
fwrite(rows,file.path(out,'within_subtype_donor_gene_expression.csv'))
fwrite(rows[,.(mean_OA=mean(logCPM[group=='OA']),mean_H=mean(logCPM[group=='Healthy']),delta=mean(logCPM[group=='OA'])-mean(logCPM[group=='Healthy'])),by=.(gene,subtype)],file.path(out,'within_subtype_gene_contrasts.csv'))

# Curves use one common H-low origin, not separate centering per simulation set.
grid<-sort(unique(c(seq(.02,.30,.01),flo,fhi)));curvefun<-function(Mac,Rest,hh=H,oo=O){mh<-rowMeans(Mac[,hh,drop=FALSE]);mo<-rowMeans(Mac[,oo,drop=FALSE]);bg<-rowMeans(Rest[,hh,drop=FALSE]);C<-sapply(c('H','OA'),function(state)sapply(grid,function(f)4000*((1-f)*bg+f*if(state=='H')mh else mo)),simplify=FALSE);C<-do.call(cbind,C);S<-score_log(logcpm(C))$scores;baseline<-score_log(logcpm(matrix(4000*((1-flo)*bg+flo*mh),ncol=1,dimnames=list(rownames(Mac),NULL))))$scores[1,];S<-sweep(S,2,baseline,'-');data.table(f=rep(grid,2),profile=rep(c('H','OA'),each=length(grid)),as.data.table(S))}
curves<-curvefun(Ms,Rs);fwrite(curves,file.path(out,'Fig5C_mixing_curves.csv'))
curve_lodo<-rbindlist(lapply(donors,function(d){z<-curvefun(Ms,Rs,setdiff(H,d),setdiff(O,d));z[,omitted_donor:=d];z}));fwrite(curve_lodo,file.path(out,'mixing_curve_LODO.csv'))

# Specified sensitivities: natural subtype mix, common-support major subtypes, equal RNA,
# alternative fixed background, old fraction anchors, and reference-anchored TMM.
sens<-list(primary_fixed8=ef,natural_subtype_mix=evaluate(Mns,Rs)$effects,major5_fixed_weights=evaluate(Mmaj,Rs)$effects,equal_RNA_per_cell=evaluate(Meq,Req)$effects,legacy_fraction_anchors=evaluate(Ms,Rs,lo=.068,hi=.185)$effects)
for(bgname in c('OA_background','balanced_all_background')){bg<-if(bgname=='OA_background')rowMeans(Rs[,O])else rowMeans(Rs);C<-cell_counts(rowMeans(Ms[,H]),rowMeans(Ms[,O]),bg);sens[[bgname]]<-effects(score_log(logcpm(C))$scores)}
fullC<-cell_counts(rowMeans(mac_std[,H]),rowMeans(mac_std[,O]),rowMeans(rest[,H]));nf<-vapply(seq_len(ncol(fullC)),function(j){dd<-DGEList(cbind(fullC[,1],fullC[,j]));dd<-normLibSizes(dd,method='TMM',refColumn=1);dd$samples$norm.factors[2]/dd$samples$norm.factors[1]},numeric(1));TE<-log2(sweep(fullC,2,colSums(fullC)*nf,'/')*1e6+1);sens$anchored_TMM<-effects(score_log(TE)$scores)
fwrite(rbindlist(sens,idcol='scenario'),file.path(out,'factorial_sensitivity_scenarios.csv'));fwrite(data.table(condition=colnames(fullC),TMM_factor=nf,library_size=colSums(fullC)),file.path(out,'TMM_sensitivity_factors.csv'))
# Fixed-background marker check: denominator-only profile effect for COL1A1.
den<-fullC;den['COL1A1',]<-4000*rep(c(1-flo,1-fhi),2)*mean(rest['COL1A1',H]);# counterfactual numerator; use the actual condition's total RNA below
rawlog<-log2(sweep(den['COL1A1',,drop=FALSE],2,colSums(fullC),'/')*1e6+1)
fwrite(data.table(condition=colnames(fullC),actual_logCPM=logcpm(fullC)['COL1A1',],denominator_only_logCPM=as.numeric(rawlog)),file.path(out,'COL1A1_fixed_background_check.csv'))

# Cross-language verification inputs: target+OTHER profiles preserve exact raw library mass.
profile_export<-rbindlist(lapply(list(standardized8=Ms,natural= Mns,major5=Mmaj,equalRNA=Meq,rest=Rs,rest_equalRNA=Req),function(A){as.data.table(A,keep.rownames='gene')}),idcol='profile_type')
fwrite(profile_export,file.path(out,'donor_count_profiles_target_plus_OTHER.csv.gz'))
fwrite(as.data.table(small(o$subtype_sum),keep.rownames='gene'),file.path(out,'source_subtype_sum_target_plus_OTHER.csv.gz'))
fwrite(as.data.table(small(o$lineage_sum),keep.rownames='gene'),file.path(out,'source_lineage_sum_target_plus_OTHER.csv.gz'))
fwrite(as.data.table(small(o$clean_mac_sum),keep.rownames='gene'),file.path(out,'source_clean_mac_sum_target_plus_OTHER.csv.gz'))
fwrite(data.table(sample=donors,n_cells=o$clean_mac_n[donors]),file.path(out,'clean_mac_counts_by_donor.csv'))
saveRDS(list(primary=primary,primary_effects=ef,effect_summary=long,subtype_weights=w,major_subtypes=major,profiles_small=profile_export,donor_scores=donexpr,curves=curves,sensitivity=sens,anchors=c(low=flo,high=fhi)),file.path(out,'Fig5_revision_objects.rds'))
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'))
cat('ANCHORS',flo,fhi,'SUBTYPE_MIN_CELLS',min(ntab),'MAJOR',major,'\n');cat('PRIMARY_EFFECTS\n');print(ef);cat('DONOR_CHECKS\n');print(pex);cat('LODO_CORE\n');print(long[variable%in%c('Composite','OLR1','PYCARD','CD68')]);cat('SENSITIVITY_CORE\n');print(rbindlist(sens,idcol='scenario')[variable%in%c('Composite','OLR1','PYCARD','CD68')]);cat('VALIDATION_OK\n')
