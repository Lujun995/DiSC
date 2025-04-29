perm.fdr.adj <- function (F0, Fp) {
  ord <- order(F0, decreasing = T)
  F0 <- F0[ord]
  perm.no <- ncol(Fp)
  Fp <- as.vector(Fp)
  Fp <- Fp[!is.na(Fp)]
  Fp <- sort(c(Fp, F0), decreasing = F)
  n <- length(Fp)
  m <- length(F0)
  FPN <- (n + 1) - match(F0, Fp) - 1:m
  p.adj.fdr <- FPN / perm.no / (1:m)
  
  # Impute 0s - pseudo-ct 0.5
  p.adj.fdr[p.adj.fdr == 0] <- 0.5 / perm.no / which(p.adj.fdr == 0)
  
  p.adj.fdr <- pmin(1, rev(cummin(rev(p.adj.fdr))))[order(ord)]
  
  return(p.adj.fdr)
}