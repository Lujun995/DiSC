DiSC <- function (data.mat, cell.ind, metadata, outcome, covariates = NULL,
                  cell.id = "cell_id", individual.id = "individual",
                  perm.no = 999, features = c('prev', 'nzm', 'nzsd'),
                  verbose = TRUE, sequencing.data = TRUE) {
  #data.mat: data matrix. Genes/variables are in rows, cells are in columns,
  # colnames are `cell_id`
  #cell.ind: a data frame linking `cell_id` to `individual_id`
  #metadata: a dataframe including an 'individual_id', an outcome of interest,
  # and covariates
  #cell.id: the variable name of cell ids in `cell.ind`
  #individual.id: the variable name of the individual id variables in `cell.ind`
  # and `metadata`
  #perm.no: number of permutations
  #features: feature extracted from the data matrix for each individual
  #verbose: do we need to print the progresses?
  #sequencing.data: is the data.mat a sequencing data matrix
  # (e.g. scRNA sequencing data/count data)?

  inds <- unique(cell.ind[[individual.id]])
  #match(target, df$name), order metadata according to inds
  metadata <- metadata[base::match(inds, metadata[[individual.id]]), ]

  cell.list <- list()
  for (ind in inds) {
    cell.list[[ind]] <- cell.ind[cell.ind[[individual.id]] == ind, ][[cell.id]]
  }

  if(sequencing.data){
    depth <- colSums(data.mat)
    data.mat <- t(t(data.mat) / depth)
    log_md_depth <- numeric(length = length(inds))
    names(log_md_depth) <- inds
    for(ind in inds)
      log_md_depth[ind] <- log(median(depth[cell.list[[ind]]]))
    metadata$log_md_depth <- log_md_depth
  } else {
    depth <- log_md_depth <- NULL
  }

  yy.list <- list()
  for (feature in features) {
    yy.list[[feature]] <- matrix(NA, nrow(data.mat), length(inds),
                                 dimnames = list(rownames(data.mat), inds))
  }
  # feature extraction
  for (ind in inds) {
    expr <- data.mat[, cell.list[[ind]]]
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
    if('nzm^1' %in% features){
      yy.list[['nzm^1']][, ind] <- rowMeans(expr)
    }
    if('nzsd^1' %in% features){
      yy.list[['nzsd^1']][, ind] <-
        apply(expr, MARGIN = 1,
              FUN = function(x){
                temp = sd(x[x > 0])
                if(is.na(temp)) return(0) else
                  return(temp)}
        )
    }
    # pretty balanced. Some genes have highest Fs in "prev", while others in
    #"nzm" or "nzsd"
  }

  obj <- basic.func(metadata, yy.list,
                    grp.name = outcome,
                    adj.name = c(covariates, ifelse(sequencing.data,
                                                    "log_md_depth", NULL)),
                    stats.combine.func = base::max,
                    perm.no = perm.no, verbose = verbose)

  return(obj)
}
