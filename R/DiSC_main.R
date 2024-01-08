DiSC <- function (ct.mat, cell.ind, metadata, outcome, covariates = NULL,
                  stats.combine.func = max, cell.id = "cell_id",
                  individual.id = "individual", perm.no = 999,
                  features = c('prev', 'nzm', 'nzsd'), verbose = TRUE, ...) {
  #ct.mat: sequencing count matrix. Genes are in rows, colnames are `cell_id`
  #cell.ind: a data frame linking `cell_id` to `individual_id`
  #metadata: a dataframe including an 'individual_id', an outcome of interest,
  # and covariates
  #cell.id: the variable name of cell ids in `cell.ind`
  #individual.id: the variable name of the individual id variables in `cell.ind` and `metadata`

  depth <- colSums(ct.mat)
  ct.mat <- t(t(ct.mat) / colSums(ct.mat))

  inds <- unique(cell.ind[[individual.id]])
  cell.list <- list()
  for (ind in inds) {
    cell.list[[ind]] <- cell.ind[cell.ind[[individual.id]] == ind, ][[cell.id]]
  }

  yy.list <- list()

  for (feature in features) {
    yy.list[[feature]] <- matrix(NA, nrow(ct.mat), length(inds),
                                 dimnames = list(rownames(ct.mat), inds))
  }

  for (ind in inds) {
    # feature extraction
    expr <- ct.mat[, cell.list[[ind]]]
    if('prev' %in% features){
      prv <- rowMeans(expr == 0)
      prv <- log(prv / (1 - prv))
      prv[prv > 7] <- 7
      prv[prv < -7] <- -7 #prv < 0.001
      yy.list[['prev']][, ind] <- prv
    }
    if('nzm' %in% features){
      yy.list[['nzm']][, ind] <- sqrt(rowMeans(expr))
    }
    if('nzsd' %in% features){
      temp <- expr
      temp[temp == 0] <- NA
      yy.list[['nzsd']][, ind] <-
        sqrt(apply(expr, MARGIN = 1,
                   FUN = function(x){
                     temp = sd(x[x > 0])
                     if(is.na(temp)) return(0) else
                       return(temp)}
      ))}
    if('sd' %in% features){
      yy.list[['sd']][, ind]   <- sqrt(rowSds(expr))
    }
    # pretty balanced. Some genes have highest Fs in "prev", while others in "nzm" or "nzsd"
  }

  #match(target, df$name)
  metadata <- metadata[base::match(inds, metadata[[individual.id]]), ]

  obj <- basic.func(metadata, yy.list,
                    grp.name = outcome,
                    adj.name = covariates,
                    stats.combine.func = stats.combine.func,
                    perm.no = perm.no, verbose = verbose,
                    ...)

  return(obj)
}
