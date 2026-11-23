Sys.setlocale('LC_ALL','English_United States.utf8')
suppressPackageStartupMessages(library(data.table))
M<-Sys.getenv('OA_REPO_ROOT',unset=getwd());S<-file.path(M,'data/components/macrophages');out<-file.path(M,'results/donor_audit');dir.create(out,recursive=TRUE,showWarnings=FALSE);A<-file.path(out,'analysis');dir.create(A,showWarnings=FALSE);dir.create(file.path(out,'figures'),showWarnings=FALSE)
writeLines(c('# 核心分解复核预设','从原始target counts独立重建，复算原始2分量和供体均值+协方差3分量。检查全部20种3对3标签分配、6次供体剔除、全部9种2对2子样本配对、6次模块剔除、44次单基因剔除、模块分项及预定义亚群合并。不同评分尺度仅检查结论适用范围，不按结果选择最佳评分。20种标签检验以供体为单位，不以细胞为单位。现有分解与论文不覆盖。'),file.path(out,'复核预设.md'))
o<-readRDS(file.path(S,'inputs/GSE216651_target_counts.rds'));md<-as.data.table(o$metadata);C<-as.matrix(o$counts);X<-log1p(sweep(C,2,o$library_total,'/')*1e4)
stopifnot(ncol(C)==7492,nrow(md)==7492,identical(colnames(C),md$cell),!anyDuplicated(md$cell),all(is.finite(X)),all(C>=0),all(o$library_total>0))
H<-sort(unique(md[group=='Healthy',sample]));O<-sort(unique(md[group=='OA',sample]));ids<-c(H,O);stopifnot(length(H)==3,length(O)==3)
gs<-fread(file.path(S,'inputs/gene_set_mitoxyperilysis.csv'));setnames(gs,1,'module');mods<-c('M2_mito_damage','M3_antiox_buffer','M4_mtorc2','M4_mtorc1','M5_cytoskeleton','M6_membrane_lysis');g<-gs[module%in%mods];stopifnot(uniqueN(g$gene)==44,all(g$gene%in%rownames(X)))
cw<-1/(length(ids)*as.numeric(table(md$sample)[md$sample]));mu<-drop(X%*%cw);sdv<-sqrt(drop((X-mu)^2%*%cw));sdv[sdv<1e-12]<-1;Z<-(X-mu)/sdv
score<-function(Z,gm=g,modules=mods){ms<-sapply(modules,function(m){q<-gm[module==m];colMeans(Z[q$gene,,drop=FALSE]*q$direction)});rowMeans(ms)}
aggregate_score<-function(y,map=NULL){
 dd<-copy(md[,.(sample,group,subtype)]);if(!is.null(map))dd[,subtype:=map[as.character(subtype)]];stopifnot(!anyNA(dd$subtype));dd[,value:=y]
 dd<-dd[,.(e=mean(value),n=.N),by=.(sample,group,subtype)];dd[,w:=n/sum(n),by=sample];dd
}
decompose<-function(dd,hh=H,oo=O){
 kk<-sort(unique(dd$subtype));part<-function(donors){rbindlist(lapply(kk,function(k){q<-dd[subtype==k];q<-q[match(donors,q$sample)];stopifnot(!anyNA(q$sample));w<-mean(q$w);e<-mean(q$e);cv<-mean((q$w-w)*(q$e-e));data.table(subtype=k,w=w,e=e,cov=cv,joint_e=mean(q$w*q$e)/w)}))}
 h<-part(hh);a<-part(oo);p<-data.table(subtype=kk,w_H=h$w,w_OA=a$w,e_H=h$e,e_OA=a$e,cov_H=h$cov,cov_OA=a$cov,joint_e_H=h$joint_e,joint_e_OA=a$joint_e)
 p[,`:=`(composition=(w_OA-w_H)*(joint_e_H+joint_e_OA)/2,within_profile=(w_H+w_OA)*(joint_e_OA-joint_e_H)/2,composition_donor_mean=(w_OA-w_H)*(e_H+e_OA)/2,expression_donor_mean=(w_H+w_OA)*(e_OA-e_H)/2,covariance_change=cov_OA-cov_H)]
 observed<-dd[,.(s=sum(w*e)),by=sample];delta<-mean(observed$s[match(oo,observed$sample)])-mean(observed$s[match(hh,observed$sample)])
 totals<-p[,lapply(.SD,sum),.SDcols=c('composition','within_profile','composition_donor_mean','expression_donor_mean','covariance_change')];totals[,total:=delta]
 stopifnot(abs(totals$composition+totals$within_profile-delta)<1e-10,abs(totals$composition_donor_mean+totals$expression_donor_mean+totals$covariance_change-delta)<1e-10)
 list(parts=p,total=totals)
}
y<-score(Z);dd<-aggregate_score(y);ans<-decompose(dd)
old<-fread(file.path(S,'analysis/decomposition_total.csv'));stopifnot(abs(old$composition-ans$total$composition)<1e-10,abs(old$within_profile-ans$total$within_profile)<1e-10)
fwrite(ans$parts,file.path(A,'primary_subtype_three_component.csv'));fwrite(ans$total,file.path(A,'primary_two_and_three_component.csv'));fwrite(dd,file.path(A,'donor_subtype_reconstructed.csv'))
don<-dd[,.(score=sum(w*e)),by=.(sample,group)];don[,library_median:=sapply(sample,function(d)median(o$library_total[md$sample==d]))];don[,n_cells:=as.numeric(table(md$sample)[sample])];fwrite(don,file.path(A,'donor_scores_and_depth.csv'))
allversions<-list(primary=y)
for(m in mods)allversions[[paste0('drop_module_',m)]]<-score(Z,modules=setdiff(mods,m))
for(m in mods)allversions[[paste0('module_',m)]]<-score(Z,modules=m)
eg<-unique(g[,.(gene,direction)]);stopifnot(!anyDuplicated(eg$gene));allversions$equal_unique_gene_weights<-colMeans(Z[eg$gene,,drop=FALSE]*eg$direction)
allversions$unscaled_log_expression<-score(X)
dm<-sapply(ids,function(d)rowMeans(X[,md$sample==d,drop=FALSE]));dmu<-rowMeans(dm);dsd<-apply(dm,1,sd);dsd[dsd<1e-12]<-1;allversions$donor_mean_gene_scaling<-score((X-dmu)/dsd)
versions<-rbindlist(lapply(names(allversions),function(name){z<-decompose(aggregate_score(allversions[[name]]))$total;z[,scenario:=name];z}));fwrite(versions,file.path(A,'score_and_module_sensitivity.csv'))
geneout<-rbindlist(lapply(unique(g$gene),function(omitted){gm<-g[g$gene!=omitted];active<-mods[mods%in%gm$module];q<-decompose(aggregate_score(score(Z,gm,active)))$total;q[,omitted_gene:=omitted];q}))
fwrite(geneout,file.path(A,'leave_one_gene_out.csv'))
lod<-rbindlist(lapply(ids,function(d){z<-decompose(dd,setdiff(H,d),setdiff(O,d))$total;z[,omitted_donor:=d];z}));fwrite(lod,file.path(A,'leave_one_donor_out.csv'))
two<-rbindlist(lapply(combn(H,2,simplify=FALSE),function(hh)rbindlist(lapply(combn(O,2,simplify=FALSE),function(oo){q<-decompose(dd,hh,oo)$total;q[,`:=`(healthy=paste(hh,collapse=';'),OA=paste(oo,collapse=';'))];q}))));fwrite(two,file.path(A,'all_nine_two_vs_two_subsets.csv'))
perms<-combn(ids,3,simplify=FALSE);null<-rbindlist(lapply(seq_along(perms),function(i){hh<-perms[[i]];q<-decompose(dd,hh,setdiff(ids,hh))$total;q[,`:=`(allocation=i,healthy=paste(hh,collapse=';'))];q}));fwrite(null,file.path(A,'all20_donor_label_permutations.csv'))
ptest<-rbindlist(lapply(setdiff(names(ans$total),'total'),function(v)data.table(term=v,observed=ans$total[[v]],exact_two_sided_p=mean(abs(null[[v]])>=abs(ans$total[[v]])-1e-12))))
ptest<-rbind(ptest,data.table(term='total',observed=ans$total$total,exact_two_sided_p=mean(abs(null$total)>=abs(ans$total$total)-1e-12)));fwrite(ptest,file.path(A,'exact_permutation_p.csv'))
subtypes<-unique(md$subtype);base<-setNames(as.character(subtypes),subtypes)
maps<-list(original8=base)
m<-base;m[grepl('Resident',names(m))]<-'Resident_merged';maps$resident_merged6<-m
m<-maps$resident_merged6;m[names(m)%in%c('Mac_OLR1','LAM_SPP1')]<-'Remodeling_merged';maps$resident_remodeling_merged5<-m
m<-base;m[names(m)%in%c('Resident_CCL2','Resident_MARCO','Resident_LYVE1')]<-'Resident';m[names(m)%in%c('Mac_OLR1','LAM_SPP1')]<-'Remodeling';m[!m%in%c('Resident','Remodeling')]<-'Other';maps$broad3<-m
merges<-rbindlist(lapply(names(maps),function(name){q<-decompose(aggregate_score(y,maps[[name]]))$total;q[,scenario:=name];q}));fwrite(merges,file.path(A,'coarsened_subtype_sensitivity.csv'));fwrite(rbindlist(lapply(names(maps),function(n)data.table(scenario=n,original_subtype=names(maps[[n]]),merged_subtype=unname(maps[[n]])))),file.path(A,'subtype_merge_dictionary.csv'))
# Within-donor depth-matched scoring: common depth strata across all six donors.
breaks<-unique(quantile(o$library_total,seq(0,1,length.out=5)));bin<-cut(o$library_total,breaks,include.lowest=TRUE,labels=FALSE)
depth<-md[,.(cell,sample,group,subtype)];depth[,`:=`(bin=bin,value=y)]
depth_summary<-depth[,.(score=mean(value),n=.N),by=.(sample,group,bin)];common<-depth_summary[,.(donors=uniqueN(sample),min_n=min(n)),by=bin][donors==6&min_n>=20,bin]
if(length(common)){
 ds<-depth_summary[bin%in%common,.(score=mean(score),bins=.N),by=.(sample,group)];fwrite(ds,file.path(A,'common_depth_strata_donor_scores.csv'))
}
fwrite(g[,.(gene,direction,weight=direction/(length(mods)*.N)),by=module],file.path(A,'effective_gene_weights.csv'))
writeLines(capture.output(sessionInfo()),file.path(A,'R_sessionInfo.txt'))
cat('PRIMARY\n');print(ans$total);cat('MODULES AND SCALING\n');print(versions);cat('GENE LOO OPPOSITE',sum(geneout$composition<0&geneout$within_profile>0),'/',nrow(geneout),'\n');cat('ALL2v2 OPPOSITE',sum(two$composition<0&two$within_profile>0),'/',nrow(two),'\n');print(ptest);print(merges)
pdf(file.path(out,'figures/Core_decomposition_audit.pdf'),width=9,height=6)
par(mar=c(5,10,3,2));vs<-versions[scenario=='primary'|grepl('^drop_module_',scenario)|scenario=='equal_unique_gene_weights'];yy<-seq_len(nrow(vs));plot(vs$composition,yy,xlim=c(min(vs$composition)-.01,max(vs$within_profile)+.13),yaxt='n',ylab='',xlab='Decomposition component',pch=19,col='#0072B2',main='Score construction sensitivity');points(vs$within_profile,yy,pch=17,col='#D55E00');axis(2,at=yy,labels=sub('drop_module_','Without ',vs$scenario),las=1,cex.axis=.75);abline(v=0,lty=2);legend('topright',c('Composition','Within-profile'),pch=c(19,17),col=c('#0072B2','#D55E00'),bty='n')
par(mar=c(5,5,3,2));bars<-as.numeric(ans$total[,.(composition_donor_mean,expression_donor_mean,covariance_change)]);barplot(bars,ylim=extendrange(bars,f=.15),names.arg=c('Composition','Expression','Donor covariance'),col=c('#0072B2','#D55E00','#009E73'),ylab='Score difference',main='Donor-mean three-component decomposition');abline(h=0)
dev.off()
