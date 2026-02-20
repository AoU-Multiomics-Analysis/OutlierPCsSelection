
ComputePCs <- function(SamplesByGenes) {
PCARes <- PCAtools::pca(samplesByGenes %>% t())
PCARes
}


ComputePCBetaValues <- function(PCScores,SamplesByGenes) {
    Z <- as.matrix(PCScores)
    ExpressionMatrix <- data.matrix(SamplesByGenes)
    stopifnot(nrow(Z) == nrow(SamplesByGenes))
    CenteredExpression <- sweep(ExpressionMatrix, 2,
                            colMeans(ExpressionMatrix, na.rm = TRUE),
                            FUN = "-")    
    Beta <- crossprod(Z,CenteredExpression)/colSums(Z^2)
    Beta   
}

ZScores_Kgrid <- function(Z, BetaPCs, SamplesByGenes, K_grid,
                          gene_block = 2000) {

  Z <- as.matrix(Z)       # n × P
  B <- as.matrix(BetaPCs) # P × G
  Y <- as.matrix(SamplesByGenes)

  stopifnot(nrow(Z) == nrow(Y))
  stopifnot(ncol(Z) == nrow(B))
  stopifnot(ncol(Y) == ncol(B))

  n <- nrow(Z)
  P <- ncol(Z)
  G <- ncol(B)

  sample_ids <- rownames(Y) %||% rownames(Z)
  gene_ids   <- colnames(Y)

  if (is.null(sample_ids)) sample_ids <- seq_len(n)
  if (is.null(gene_ids))   gene_ids   <- seq_len(G)

  K_grid <- sort(unique(K_grid))
  stopifnot(all(K_grid >= 0), all(K_grid <= P))
  Kmax <- max(K_grid)

  # Center expression
  Y0 <- sweep(Y, 2, colMeans(Y, na.rm = TRUE), "-")

  # Preallocate output with dimnames
  z_list <- lapply(K_grid, function(k) {
    matrix(NA_real_,
           nrow = n, ncol = G,
           dimnames = list(sample_ids, gene_ids))
  })
  names(z_list) <- as.character(K_grid)

  for (g0 in seq(1, G, by = gene_block)) {

    g1 <- min(G, g0 + gene_block - 1)

    Bb <- B[, g0:g1, drop = FALSE]        # P × Gb
    Yb <- Y0[, g0:g1, drop = FALSE]       # n × Gb

    fit_block <- matrix(0, nrow = n, ncol = ncol(Bb))

    # K = 0 case
    if (0 %in% K_grid) {
      resid <- Yb
      df <- n - 1
      sigma <- sqrt(colSums(resid^2) / df)
      z_list[["0"]][, g0:g1] <- sweep(resid, 2, sigma, "/")
    }

    for (k in 1:Kmax) {

      fit_block <- fit_block +
        Z[, k, drop = FALSE] %*% t(Bb[k, , drop = TRUE])

      if (k %in% K_grid) {

        resid <- Yb - fit_block
        df <- n - 1 - k
        if (df <= 1) stop("Not enough df left; reduce K.")

        sigma <- sqrt(colSums(resid^2) / df)

        z_list[[as.character(k)]][, g0:g1] <-
          sweep(resid, 2, sigma, "/")
      }
    }
  }

  z_list
}

