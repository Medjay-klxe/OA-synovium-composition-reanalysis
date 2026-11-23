invisible(Sys.setlocale('LC_ALL','Chinese (Simplified)_China.utf8'))

args <- commandArgs(TRUE)
stopifnot(length(args)==4)
old_script <- args[1]; output <- args[2]; batch <- as.integer(args[3]); B <- as.integer(args[4])
prefix <- readLines(old_script, encoding='UTF-8', warn=FALSE)
prefix <- sub('; baseout<-.*$', '', prefix)
prefix <- prefix[seq_len(grep('^# Simulated reference discrepancy',prefix)[1]-1)]
eval(parse(text=prefix[!startsWith(prefix,"a<-commandArgs")]))
suppressPackageStartupMessages(library(GSVA))
set.seed(20261130L+batch)
donors <- p$donor_metadata$sample
stopifnot(length(donors)==6L,all(donors%in%md$sample))
build_heldout <- function(omit) {
  pools <- lapply(cls,function(cell_type) {
    if(cell_type=='Macrophage') return(rowSums(p$clean_mac_sum[,setdiff(donors,omit),drop=FALSE]))
    rowSums(p$lineage_sum[,md[!sample %in% omit & cl==cell_type,profile_id],drop=FALSE])
  })
  m <- do.call(cbind,pools)
  rownames(m) <- rownames(p$lineage_sum); colnames(m) <- cls
  stopifnot(all(colSums(m)>0))
  sweep(m,2,colSums(m),'/')[mk,,drop=FALSE]*1e6
}
references <- c(list(all6=ref[mk,,drop=FALSE]),setNames(lapply(donors,build_heldout),paste0('omit_',donors)))
estimate <- function(raw,mb,observed) {
  avg <- rowMeans(raw)
  calibrated <- raw*ifelse(avg>1,mb/avg,1)
  stopifnot(all(colSums(calibrated)>0))
  Sfit <- sweep(calibrated,2,colSums(calibrated),'/')
  scale_gene <- apply(Sfit,1,sd); scale_gene[scale_gene<1e-8] <- 1
  center <- rowMeans(Sfit)
  A <- (Sfit-center)/scale_gene
  D <- crossprod(A)+diag(1e-8,ncol(A))
  d <- crossprod(A,(observed-center)/scale_gene)
  constraints <- cbind(rep(1,ncol(A)),diag(ncol(A)))
  vapply(seq_len(ncol(d)),function(j) {
    sol <- solve.QP(D,d[,j],constraints,c(1,rep(0,ncol(A))),meq=1)$solution
    stopifnot(abs(sum(sol)-1)<1e-6,min(sol)>-1e-6)
    sol[match('Macrophage',cls)]
  },numeric(1))
}
score <- function(C) {
  d <- DGEList(C)
  d <- d[filterByExpr(d,group=factor(label)),,keep.lib.sizes=FALSE]
  E <- cpm(normLibSizes(d,method='TMM'),log=TRUE,prior.count=1)
  use <- lapply(sets,intersect,y=rownames(E))
  stopifnot(length(use)==39L,min(lengths(use))>=3L)
  Z <- t(scale(t(E))); Z[!is.finite(Z)] <- 0
  ssg <- t(gsva(ssgseaParam(E,use,minSize=3,alpha=.25,normalize=TRUE,verbose=FALSE),verbose=FALSE))
  stopifnot(setequal(colnames(ssg),names(sets)))
  list(z_mean=scale(sapply(use,function(g) colMeans(Z[g,,drop=FALSE]))), ssGSEA=scale(ssg[,names(sets),drop=FALSE]))
}
fit <- function(Y,covariate=NULL) {
  X <- cbind(Intercept=1,label=label,covariate)
  model <- lm.fit(X,Y)
  stopifnot(model$rank==ncol(X),all(is.finite(model$coefficients)))
  beta <- model$coefficients[2,]
  se <- sqrt(colSums(model$residuals^2)/model$df.residual*solve(crossprod(X))[2,2])
  list(beta=beta,p=2*pt(-abs(beta/se),model$df.residual))
}
rem <- setdiff(colnames(bulkcp),meta$sample)
bg <- rowSums(bulkcp[mk,rem,drop=FALSE])
out <- list(); qc <- list()
scenarios <- c('no_change','composition_only','expression_only','both')
cat('READY',batch,'GSVA',as.character(packageVersion('GSVA')),'\n')
for(rep in seq_len(B)) {
  for(sc in scenarios) {
    w <- if(sc %in% c('no_change','expression_only')) prop[sample(nrow(prop)),,drop=FALSE] else prop
    has_expr <- sc %in% c('expression_only','both')
    Y <- G%*%t(w)
    if(has_expr) Y <- Y+tcrossprod(GO[,'Macrophage']-G[,'Macrophage'],w[,'Macrophage']*label)
    mu <- sweep(Y,2,L,'*')
    C <- matrix(rnbinom(length(mu),mu=as.vector(mu),size=50),nrow=nrow(mu),dimnames=list(rownames(G),meta$sample))
    stopifnot(all(is.finite(C)),all(colSums(C[mk,])>0))
    mb <- (bg+rowSums(sweep(C[mk,,drop=FALSE],2,colSums(C),'/')*1e6))/(length(rem)+ncol(C))
    obs <- sweep(C[mk,,drop=FALSE],2,colSums(C[mk,,drop=FALSE]),'/')
    delta_mass <- has_expr*label*(sum(GO[mk,'Macrophage'])-sum(G[mk,'Macrophage']))
    truth <- w[,'Macrophage']*(sum(G[mk,'Macrophage'])+delta_mass)/(as.numeric(w%*%colSums(G[mk,,drop=FALSE]))+w[,'Macrophage']*delta_mass)
    fitted_prop <- lapply(references,estimate,mb=mb,observed=obs)
    y <- score(C)
    for(ref_id in names(fitted_prop)) {
      v <- fitted_prop[[ref_id]]
      qc[[length(qc)+1L]] <- data.table(batch,rep,scenario=sc,reference=ref_id,rmse=sqrt(mean((v-truth)^2)),r=cor(v,truth),shift_from_all6=sqrt(mean((v-fitted_prop$all6)^2)))
    }
    for(method in names(y)) {
      crude <- fit(y[[method]])
      for(ref_id in c(names(fitted_prop),'oracle')) {
        adjfit <- fit(y[[method]],if(ref_id=='oracle')truth else fitted_prop[[ref_id]])
        out[[length(out)+1L]] <- data.table(batch,rep,scenario=sc,reference=ref_id,method,gene_set=names(sets),beta_before=crude$beta,beta_after=adjfit$beta,p_before=crude$p,p_after=adjfit$p,q_before=p.adjust(crude$p,'BH'),q_after=p.adjust(adjfit$p,'BH'))
      }
    }
  }
  if(rep %% 5L==0L || rep==B) {cat('DONE',batch,rep,'/',B,'\n');flush.console()}
}
fwrite(rbindlist(out),file.path(output,sprintf('fits_%d.csv.gz',batch)))
fwrite(rbindlist(qc),file.path(output,sprintf('reference_%d.csv',batch)))
