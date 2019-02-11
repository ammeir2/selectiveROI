library(selectiveROI)
args <- commandArgs(TRUE)
eval(parse(text=args[[1]]))
seed <- as.numeric(seed)


runSim <- function(config) {
  p <- 100
  musd <- config[["musd"]]
  rho <- config[["rho"]]
  th <- config[["threshold"]]
  reps <- config[["reps"]]

  simResults <- vector(reps, mode = "list")
  times <- matrix(nrow = reps, ncol = 2)
  colnames(times) <- c("gibbsTime", "sliceTime")
  errors <- matrix(nrow = reps, ncol = 2)
  colnames(errors) <- c("gibbs", "slice")
  for(m in 1:reps) {
    sigma <- rho^as.matrix(dist(1:p))
    mu <- predict(smooth.spline(rnorm(p, sd = musd), df = p/5))$y
    y <- as.numeric(mvtnorm::rmvnorm(1, mean = mu, sigma = sigma))
    threshold <- th / sqrt(diag(sigma))
    select <- abs(y) > threshold
    if(!any(select)) {
      m <- m - 1
      next
    }

    clusters <- list()
    active <- FALSE
    for(i in 1:length(y)) {
      if(!active) {
        if(select[i]) {
          start <- i
          active <- TRUE
        }
      }
      if(active) {
        if(!select[i] | sign(y[i]) != sign(y[start])) {
          end <- i - 1
          active <- FALSE
          clusters[[length(clusters) + 1]] <- start:end
        }
      }
    }

    # Keeping only largest cluster ------
    clustsize <- sapply(clusters, length)
    keep <- clusters[[which.max(clustsize)[1]]]
    keep <- c(min(keep) - 1, keep, max(keep) + 1)
    keep <- keep[keep >= 1 & keep <= p]
    y <- y[keep]
    mu <- mu[keep]
    sigma <- sigma[keep, keep]
    select <- select[keep]
    threshold <- threshold[keep]

    newmle <- NULL
    try(sliceTime <- system.time(newmle <- roiMLE(y, sigma, threshold,
                                 coordinates = matrix(1:length(y), ncol = 1),
                                 selected = select,
                                 stepSizeCoef = 0.15, stepRate = 0.65,
                                 delay = 100, maxiter = 4000,
                                 sampPerIter = 100,
                                 assumeConvergence = 2000,
                                 imputeBoundary = "neighbors"))[[3]])
    if(is.null(newmle)) {
      m <- m - 1
      next
    }

    gibbsTime <- system.time(oldmle <- optimizeSelected(y, sigma, threshold,
                                           coordinates = matrix(1:length(y), ncol = 1),
                                           selected = select,
                                           stepSizeCoef = 0.15, stepRate = 0.65,
                                           delay = 100, maxiter = 4050,
                                           assumeConvergence = 4000,
                                           probMethod = "all",
                                           imputeBoundary = "neighbors"))[[3]]
    result <- list()
    result$params <- c(musd = musd, rho = rho, threshold = th)
    result$mu <- mu
    result$selected <- select
    result$times <- c(gibbs = gibbsTime, slice = sliceTime)
    result$observed <- y
    result$newmle <- newmle$conditional
    result$oldmle <- oldmle$conditional

    # Printing size
    print(c(m = m, size = sum(select)))

    # Printing times
    times[m, ] <- result$times
    print(colMeans(times[1:m, , drop = FALSE]))

    # Printing bias
    mu <- mu * sign(y[2])
    slice <- newmle$conditional * sign(y[2])
    gibbs <- oldmle$conditional * sign(y[2])
    errors[m, ] <- c(gibbs = mean(gibbs[select]) - mean(mu[select]),
                      slice = mean(slice[select]) - mean(mu[select]))
    print(colMeans(errors[1:m, , drop = FALSE]))

    # Saving results
    simResults[[m]] <- result
  }

  return(simResults)
}

params <- expand.grid(p = 100,
                      musd = c(1, 2, 4),
                      rho = c(0, 0.9),
                      threshold = c(1.5, 3),
                      reps = 5)
set.seed(seed)
result <- apply(params, 1, runSim)
filename <- paste("results/sliceComparison_C_", seed, ".rds", sep = "")
saveRDS(result, file = filename)


