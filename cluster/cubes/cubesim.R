run.sim <- function(config, noise_type ="sim", noise_dat = NULL) {
  snr <- config[["snr"]]
  spread <- config[["spread"]]
  BHlevel <- config[["BHlevel"]]
  replications <- config[["replications"]]
  slack <- config[["slack"]]

  if(noise_type == "sim" ) {
    rho <- config[["rho"]]
    grp_size <- -1
  }

  I <- 11
  J <- 11
  K <- 9
  dims <- c(I, J, K)

  coordinates <- expand.grid(i = 1:I, j = 1:J, k = 1:K)

  if (noise_type == "sim") {
    covariance <- rho^as.matrix(dist(coordinates[, 1:3], method = "euclidean",
                                     diag = TRUE, upper = TRUE))
    covEigen <- eigen(covariance)
    sqrtCov <- covEigen$vectors %*% diag(sqrt(covEigen$values)) %*% t(covEigen$vectors)
  }

  simresults <- list()
  simcover <- matrix(nrow = replications, ncol = 2)
  colnames(simcover) <- c("cover", "power")
  origCov <- covariance
  for(rep in 1:replications) {
    selected <- FALSE
    maxsize <- 0
    covariance <- origCov
    while(is.na(maxsize) | maxsize < 2) {
      if (noise_type == "sim") {
        coordinates <- generateArrayData3D(dims, sqrtCov, snr, spread)
        sds <- 1
      }

      coordinates$observed <- coordinates$observed / sds
      coordinates$signal <- coordinates$signal / sds
      coordinates$zval <- coordinates$observed
      coordinates$pval <- 2 * pnorm(-abs(coordinates$zval))
      coordinates$qval <- p.adjust(coordinates$pval, method = "bonferroni")
      coordinates$qval <- coordinates$pval  # FOR NO CORRECTION!!!!
      coordinates$selected <- coordinates$qval < BHlevel
      selected <- coordinates$selected
      if(sum(selected) >= 2) {
        threshold <- abs(qnorm(BHlevel / 2))
        coordinates$selected <- coordinates$zval > threshold
        pclusters <- findClusters(coordinates)

        coordinates$selected <- coordinates$zval < -threshold
        nclusters <- findClusters(coordinates)
        clusters <- c(pclusters, nclusters)
        coordinates$selected <- coordinates$qval < BHlevel

        maxsize <- max(sapply(clusters, function(x) sum(x$selected)))
      }
    }

    coordinates$observed <- coordinates$observed * sds
    coordinates$signal <- coordinates$signal * sds
    threshold <- abs(qnorm(BHlevel / 2)) * sds

    #covariance <- cov2cor(covariance)

    print(c(rep = rep, rho = rho, grp_size = grp_size, BHlevel = BHlevel, snr = snr, spread = spread))
    print(c(nselected = sum(selected)))
    sizes <- sapply(clusters, nrow)

    results <- list()
    iterPower <- 0
    iterCover <- 0
    weights <- 0
    slot <- 1
    mse <- 0
    msecount <- 0
    for(m in 1:length(clusters)) {
      cluster <- clusters[[m]]
      cluster <- subset(cluster, !is.na(cluster$selected))
      selected <- coordinates$selected[cluster$row]
      observed <- coordinates$observed[cluster$row]
      signal <- coordinates$signal[cluster$row]
      subCov <- covariance[cluster$row, cluster$row, drop = FALSE]

      print(c(round(m / length(clusters), 2), nrow(cluster)))

      try(mle <- roiMLE(observed, subCov, threshold[cluster$row],
                        selected = selected,
                        projected = NULL,
                        stepRate = 0.6, stepSizeCoef = 2,
                        coordinates = cluster[, 1:3],
                        tykohonovSlack = slack,
                        delay = 10,
                        assumeConvergence = 1500,
                        sampPerIter = 100,
                        maxiter = 2000,
                        nsamp = 500,
                        imputeBoundary = "neighbors"))

      w <- rep(1, sum(selected))
      w <- w / sum(w)
      conditional <- weighted.mean(mle$conditional[selected], w)
      naive <- weighted.mean(observed[selected], w)
      true <- weighted.mean(signal[selected], w)
      mse <- mse * msecount / (msecount + 1) + c(naive = (naive - true)^2, conditional = (conditional - true)^2) / (msecount + 1)
      msecount <- msecount + 1
      for(methodind in 2) {
        if(methodind == 1) {
          method <- "onesided"
        } else if(methodind == 2) {
          method <- "selected"
        }
        results[[slot]] <- list()

        # Computing the p-value based on samples from the null
        nullfit <- NULL
        try(nullfit <- roiMLE(observed, subCov, threshold[cluster$row],
                                        projected = 0,
                                        selected = selected,
                                        coordinates = cluster[, 1:3],
                                        tykohonovSlack = 0.0001,
                                        stepSizeCoef = 0,
                                        delay = 1,
                                        assumeConvergence = 1,
                                        maxiter = 2,
                                        nsamp = 5000,
                                        init = rep(0, length(observed)),
                                        imputeBoundary = "neighbors"))
        if(is.null(nullfit)) next
        naive <- mean(observed[selected])
        samp <- nullfit$sample
        nullmeans <- NULL
        try(nullmeans <- as.numeric(samp[, selected, drop = FALSE] %*% w))
        if(is.null(nullmeans)) next
        # p-value based on one sided test
        if(naive > 0) {
          pvalue <- mean(naive < nullmeans)
        } else {
          pvalue <- mean(naive > nullmeans)
        }

        # Projecting to the truth, rejecting the test here implies
        # that the CI doesn't cover the truth.
        try(profile <- roiMLE(observed, subCov, threshold[cluster$row],
                          selected = selected,
                          projected = weighted.mean(signal[selected], w),
                          stepRate = 0.65, stepSizeCoef = 0.25,
                          coordinates = cluster[, 1:3],
                          tykohonovSlack = slack,
                          delay = 50,
                          assumeConvergence = 750,
                          sampPerIter = 100,
                          maxiter = 1500,
                          nsamp = 5000,
                          imputeBoundary = "neighbors"))

        samp <- profile$sample
        profMeans <- NULL
        try(profMeans <- as.numeric(samp[, selected, drop = FALSE] %*% w))
        if(is.null(profMeans)) next
        # Two sided CI
        profPval <- 2 * min(mean(naive < profMeans), mean(naive > profMeans))

        true <- weighted.mean(signal[selected], w)
        profResult <- c(true = mean(signal[selected]), profPval = profPval, pvalue = pvalue)
        results[[slot]][[1]] <- c(snr = snr, spread = spread, method = methodind,
                                  rho = rho, BHlevel = BHlevel, grp_size = grp_size,
                                  size = sum(selected), true = true, slack = slack)
        results[[slot]][[2]] <- profResult
        results[[slot]][[3]] <- c(conditional = conditional, naive = naive, true = true)

        print(profResult)
        print(mse)

        weight <- 1
        iterPower <- iterPower + weight * (pvalue < 0.05)
        iterCover <- iterCover + weight * (profPval > 0.05)
        weights <- weights + weight

        # plot(density((rowMeans(profile$sample[, selected, drop = FALSE]))), xlim = c(-4.5, 4.5))
        # abline(v = mean(profile$conditional[selected]))
        # abline(v = mean(result$conditional[selected]), col = "blue")
        # abline(v = naive, col = "red")
        # lines(density((rowMeans(result$sample[, selected]))), col = "blue")
        # abline(v = c(threshold, - threshold), col = "grey")
        slot <- slot + 1
      }
    }

    simcover[rep, 1] <- sum(iterCover) / sum(weights)
    simcover[rep, 2] <- sum(iterPower) / sum(weights)
    print(c(cover = iterCover, power = iterPower) / weights)
    print(colMeans(simcover[1:rep, , drop = FALSE]))
    simresults[[rep]] <- results
  }

  return(simresults)
}

configurations <- expand.grid(snr = c(4, 2, 1, 0),
                              spread = c(2),
                              rho = c(0.45, 0.9),
                              BHlevel = c(0.001, 0.01),
                              replications = 10,
                              slack = 1)

simResults <- apply(configurations, 1, run.sim, noise_type ="sim")
filename <- paste("results/cubeROI_A", seed, ".rds", sep = "")
saveRDS(simResults, file = filename)
