
invisible(Sys.setlocale('LC_ALL','Chinese (Simplified)_China.utf8'))
suppressPackageStartupMessages(library(data.table))
path <- commandArgs(TRUE)[1]
load_batches <- function(prefix,ext) rbindlist(lapply(1:4,function(i) {
  f <- file.path(path,paste0(prefix,i,ext))
  if (ext=='.csv.gz') fread(text=paste(readLines(gzfile(f)),collapse='\n')) else fread(f)
}))
f <- load_batches('fits_','.csv.gz')
a <- load_batches('reference_','.csv')
nk <- c('batch','rep','scenario','reference','method')
stopifnot(nrow(f)==400L*8L*2L*39L,nrow(a)==400L*7L,
          all(is.finite(as.matrix(f[,.(beta_before,beta_after,p_before,p_after)]))),
          !anyDuplicated(f[,c(nk,'gene_set'),with=FALSE]))
check <- f[,.(count=.N,error=max(abs(q_before-p.adjust(p_before,'BH')),abs(q_after-p.adjust(p_after,'BH')))),by=nk]
stopifnot(all(check$count==39),max(check$error)<1e-12)
paired <- f[,.(n_before=uniqueN(beta_before)),by=.(batch,rep,scenario,method,gene_set)]
stopifnot(all(paired$n_before==1L))
runs <- f[,.(median_delta=median(beta_before-beta_after),sig_before=sum(q_before<.05),sig_after=sum(q_after<.05)),by=nk]
full <- runs[reference=='all6',.(batch,rep,scenario,method,base_delta=median_delta,base_sig_after=sig_after)]
paired_runs <- merge(runs,full,by=c('batch','rep','scenario','method'))
paired_runs[,delta_shift:=median_delta-base_delta]
paired_runs[,sig_shift:=sig_after-base_sig_after]
summary <- paired_runs[,.(replicates=.N,median_attenuation=median(median_delta),p025=quantile(median_delta,.025),p975=quantile(median_delta,.975),mean_sig_before=mean(sig_before),mean_sig_after=mean(sig_after),pct_none_after=100*mean(sig_after==0),paired_median_shift=median(delta_shift),paired_p95_abs_shift=quantile(abs(delta_shift),.95),pct_changed_sig_count=100*mean(sig_shift!=0)),by=.(scenario,method,reference)]
accuracy <- a[,.(median_rmse=median(rmse),median_cor=median(r),median_shift_from_full=median(shift_from_all6)),by=.(scenario,reference)]
fwrite(summary,file.path(path,'sensitivity_summary.csv'))
fwrite(accuracy,file.path(path,'composition_accuracy.csv'))
fwrite(paired_runs,file.path(path,'paired_replicate_summary.csv'))
print(summary[reference=='all6'],nrows=8)
print(summary[startsWith(reference,'omit_'),.(max_median_shift=max(abs(paired_median_shift)),max_p95=max(paired_p95_abs_shift),min_none=min(pct_none_after),max_none=max(pct_none_after)),by=.(scenario,method)])
print(accuracy[,.(rmse_full=median_rmse[reference=='all6'],rmse_min=min(median_rmse),rmse_max=max(median_rmse)),by=scenario])
cat('Verified',nrow(f),'gene set model rows;',nrow(a),'composition fits\n')
