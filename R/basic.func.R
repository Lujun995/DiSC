basic.func <-
  function (meta.dat,
            feature.list,
            grp.name,
            adj.name = NULL,
            stats.combine.func = mean,
            perm.no = 999,
            strata = NULL,
            is.fwer = FALSE,
            verbose = TRUE)
  {
    this.call <- match.call()
    
    
    feature.no <- length(feature.list)
    
    sample.no <- ncol(feature.list[[1]])
    otu.no <- nrow(feature.list[[1]])
    row.names <- rownames(feature.list[[1]])
    
    if (verbose)
      cat("The data has ",
          sample.no,
          " samples and ",
          otu.no,
          " genes will be tested!\n")
    
    
    if (!is.null(strata)) {
      strata <- factor(strata)
    }
    if (is.null(adj.name)) {
      M0 <- model.matrix( ~ 1, meta.dat)
    }
    else {
      data0 <- meta.dat[, c(adj.name), drop = FALSE]
      if (sum(is.na(data0)) != 0) {
        stop("Please remove or impute NAs in the variables to be adjusted!\n")
      }
      M0 <- model.matrix( ~ ., data0)
    }
    data1 <- meta.dat[, c(grp.name), drop = FALSE]
    if (sum(is.na(data1)) != 0) {
      stop("Please remove or impute NAs in the variable of interest!\n")
    }
    M1 <- model.matrix( ~ ., data1)[,-1, drop = FALSE]
    M01 <- cbind(M0, M1)
    qrX0 <- qr(M0, tol = 1e-07)
    Q0 <- qr.Q(qrX0)
    Q0 <- Q0[, 1:qrX0$rank, drop = FALSE]
    H0 <- (Q0 %*% t(Q0))
    qrX1 <- qr(M1, tol = 1e-07)
    Q1 <- qr.Q(qrX1)
    Q1 <- Q1[, 1:qrX1$rank, drop = FALSE]
    qrX01 <- qr(M01, tol = 1e-07)
    Q01 <- qr.Q(qrX01)
    Q01 <- Q01[, 1:qrX01$rank, drop = FALSE]
    R0 <- as.matrix(resid(lm(Q1 ~ Q0 - 1)))
    pX0 <- ncol(Q0)
    pX1 <- ncol(Q1)
    pX01 <- ncol(Q01)
    df.model <- pX01 - pX0
    df.residual <- sample.no - pX01
    

    if (verbose)
      cat("Permutation testing ...\n")

    Y <- matrix(NA, sample.no, feature.no * otu.no)
    
    for (j in 1:length(feature.list)) {
      Y[, feature.no * (0:(otu.no - 1)) + j] <- t(feature.list[[j]])
    }
    
    Y <- t(Y)
    TSS <- rowSums(Y ^ 2)
    MSS01 <- rowSums((Y %*% Q01) ^ 2)
    MSS0 <- rowSums((Y %*% Q0) ^ 2)
    MSS <- (MSS01 - MSS0)
    RSS <- (TSS - MSS01)
    getPermuteMatrix <- getFromNamespace("getPermuteMatrix",
                                         "vegan")
    perm.ind <-
      getPermuteMatrix(perm.no, sample.no, strata = strata)
    perm.no <- nrow(perm.ind)
    MRSSp <- sapply(1:perm.no, function(ii) {
      if (verbose) {
        if (ii %% 10 == 0)
          cat(".")
      }
      Rp <- R0[perm.ind[ii,], , drop = FALSE]
      Rp <- Rp - H0 %*% Rp
      qrRp <- qr(Rp, tol = 1e-07)
      Q1p <- qr.Q(qrRp)
      Q1p <- Q1p[, 1:qrRp$rank, drop = FALSE]
      MSS01p <- MSS0 + rowSums((Y %*% Q1p) ^ 2)
      MSSp <- (MSS01p - MSS0)
      RSSp <- (TSS - MSS01p)
      c(MSSp, RSSp)
    })

    unit <- feature.no * otu.no
    MSSp <- MRSSp[1:unit,]
    RSSp <- MRSSp[(unit + 1):(2 * unit),]
    RSS.m <- t(array(RSS, c(feature.no, otu.no)))

    F0.m <- t(array((MSS / df.model) / (RSS / df.residual), c(feature.no,
                                                            otu.no)))

    R2.m <- t(array(MSS / TSS, c(feature.no, otu.no)))

    F0 <- (MSS / df.model) / (RSS / df.residual)
    Fp <- (MSSp / df.model) / (RSSp / df.residual)
    F0 <- array(F0, c(feature.no, otu.no))
    Fp <- array(Fp, c(feature.no, otu.no, perm.no))
    F0 <- apply(F0, 2, stats.combine.func)
    Fp <- apply(Fp, c(2, 3), stats.combine.func)

    if (verbose)
      cat("\n")
    if (mean(is.na(F0)) >= 0.1) {
      warning("More than 10% observed F stats have NA! Please check! \n")
    }
    if (mean(is.na(Fp)) >= 0.1) {
      warning("More than 10% permuted F stats have NA! Please check! \n")
    }
    na.ind <- is.na(F0)
    F0 <- F0[!na.ind]
    Fp <- Fp[!na.ind,]
    which.nan.ind <- which(!na.ind)
    p.raw <- rowMeans(cbind(Fp, F0) >= F0)
    
    p.adj.fdr <- perm.fdr.adj(F0, Fp)
    p.raw <- na.pad(p.raw, na.ind)
    p.adj.fdr <- na.pad(p.adj.fdr, na.ind)

    names(p.raw) <- names(p.adj.fdr) <- rownames(R2.m) <- rownames(RSS.m) <- rownames(F0.m) <- row.names
    colnames(R2.m) <- colnames(F0.m) <- colnames(RSS.m) <- paste0("F", 1:feature.no)
    
    if (is.fwer) {
      p.adj.fwer <- perm.fwer.adj(F0, Fp)
      p.adj.fwer <- na.pad(p.adj.fwer, na.ind)
      names(p.adj.fwer) <- row.names
    }
    else {
      p.adj.fwer <- NULL
    }
    cat('*')
    coef.list <- NULL
    for (j in 1:feature.no) {
      coef.list[[j]] <- solve(t(M01) %*% M01) %*% t(M01) %*%
        t(feature.list[[j]])
    }
    
    if (verbose)
      cat("Completed!\n")
    return(
      list(
        call = this.call,
        R2 = R2.m,
        F0 = F0.m,
        RSS = RSS.m,
        df.model = df.model,
        df.residual = df.residual,
        coef.list = coef.list,
        p.raw = p.raw,
        p.adj.fdr = p.adj.fdr,
        p.adj.fwer = p.adj.fwer
      )
    )
  }