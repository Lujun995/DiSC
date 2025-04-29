perm.fwer.adj <- function (F0, Fp) {
  ord <- order(F0, decreasing = T)
  m <- length(F0)
  F0 <- F0[ord]
  Fp <- Fp[ord, , drop=FALSE]
  col.max <- Fp[m, ]
  p.adj.fwer <- sapply(m:1, function(i) {
    x <- F0[i]
    y <- Fp[i, ]
    col.max <<- ifelse(y > col.max, y, col.max)
    # Impute 0s
    mean(c(1, col.max >= x))
  })
  p.adj.fwer <- rev(p.adj.fwer)
  
  p.adj.fwer <- pmin(1, rev(cummin(rev(p.adj.fwer))))[order(ord)]
  return(p.adj.fwer)
}