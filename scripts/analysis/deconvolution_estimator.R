# Gene-standardized nonnegative quadratic program used by the original composition analysis.
# Sm: calibrated, marker-normalized genes x cell-compartment reference.
# y: bulk marker-normalized vector on the same genes and scale.
# Frozen composition estimates are supplied in data/components/bulk/inputs.
suppressPackageStartupMessages(library(quadprog))
fit_std <- function(Sm, y) {
  mu <- rowMeans(Sm); sdv <- apply(Sm, 1, stats::sd); sdv[sdv < 1e-8] <- 1
  Sz <- (Sm - mu) / sdv; yz <- (y - mu) / sdv
  D <- crossprod(Sz) + diag(1e-8, ncol(Sm))
  d <- as.numeric(crossprod(Sz, yz))
  k <- ncol(Sm)
  A <- cbind(rep(1, k), diag(k))
  r <- tryCatch(solve.QP(D, d, A, c(1, rep(0, k)), meq = 1), error = function(e) NULL)
  if (is.null(r)) return(rep(NA_real_, k))
  pmax(r$solution, 0)
}
