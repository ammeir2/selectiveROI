# Construct second differences Tykohonov Regularization matrix
computeTykohonov <- function(selected, coordinates) {
  coordinates <- coordinates[selected, , drop = FALSE]
  if(nrow(coordinates) == 1) return(matrix(0))

  distances <- as.matrix(dist(coordinates), method = "euclidean")
  diag(distances) <- Inf

  dims <- ncol(coordinates)
  p <- sum(selected)

  firstDiff <- matrix(0, nrow = p * 3 * 2, ncol = 3)
  secondDiff <- matrix(0, nrow = p * 3 * 3, ncol = 3)
  fsparseRow <- 1
  frow <- 1
  ssparseRow <- 1
  srow <- 1
  for(i in 1:nrow(coordinates)) {
    candidates <- distances[i, ] == 1
    scandidates <- distances[i, ] == 2
    for(j in 1:dims) {
      neighbor <- (1:nrow(coordinates))[candidates][which(coordinates[i, j] - coordinates[candidates, j] == - 1)]
      if(length(neighbor) == 0) {
        next
      }
      firstDiff[fsparseRow, ] <- c(frow, neighbor, -1)
      fsparseRow <- fsparseRow + 1
      firstDiff[fsparseRow, ] <- c(frow, i, 1)
      fsparseRow <- fsparseRow + 1
      frow <- frow + 1

      if(length(scandidates) == 0) next
      sneighbor <- (1:nrow(coordinates))[scandidates][which(coordinates[i, j] - coordinates[scandidates, j] == - 2)]
      if(length(sneighbor) == 0) {
        next
      }

      secondDiff[ssparseRow, ] <- c(srow, sneighbor, 1)
      ssparseRow <- ssparseRow + 1
      secondDiff[ssparseRow, ] <- c(srow, neighbor, -2)
      ssparseRow <- ssparseRow + 1
      secondDiff[ssparseRow, ] <- c(srow, i, 1)
      ssparseRow <- ssparseRow + 1
      srow <- srow + 1
    }
  }

  firstDiff <- firstDiff[firstDiff[, 1] != 0, ]
  firstDiff <- Matrix::sparseMatrix(i = firstDiff[, 1], j = firstDiff[, 2], x = firstDiff[, 3],
                                    dims = c(nrow(firstDiff), sum(selected)))
  secondDiff <- secondDiff[secondDiff[, 2] != 0, ]
  if(length(secondDiff) == 0) {
    secondDiff <- sparseMatrix(i = p, j = p, x = 0)
  } else {
    if(max(secondDiff[, 1]) < p | max(secondDiff[, 2] < p)) {
      secondDiff <- rbind(secondDiff, c(p, p, 0))
    }
    secondDiff <- sparseMatrix(i = secondDiff[, 1], j = secondDiff[, 2], x = secondDiff[, 3])
  }
  return(list(firstDiff = firstDiff, secondDiff = secondDiff))
}

# Adjusts tykohonov parametere
adjustTykohonov <- function(obsDiff, obsmean, mu, selected,
                            firstDiff, secondDiff,
                            tykohonvSlack, tykohonovParam) {
  mu <- mu[selected]
  meanmu <- mean(mu)
  for(i in 1:2) {
    if(tykohonovParam[i] == 0) next

    if(i == 1) {
      muDiff <- crossprod(mu, crossprod(firstDiff, mu))
    } else {
      muDiff <- crossprod(mu, crossprod(secondDiff, mu))
    }

    # if(i == 1 & any(sign(mu[1]) != sign(mu)) & mean(mu) > 0.5) {
    #   tykohonovParam[i] <- tykohonovParam[i] * 1.1
    #   next
    # }

    ratio <- muDiff / (obsDiff[i] * (tykohonvSlack))
    #if(i == 2) print(c(ratio, muDiff, obsDiff[i]))
    if(is.na(ratio)) {
      tykohonovParam[i] <- Inf
    } else if(ratio > 2) {
      tykohonovParam[i] <- tykohonovParam[i] * 1.3
    } else if(ratio > 1) {
      tykohonovParam[i] <- tykohonovParam[i] * 1.05
    } else if(ratio < 0.5) {
      tykohonovParam[i] <- tykohonovParam[i] * 0.7
    } else if(ratio < 1) {
      tykohonovParam[i] * 0.95
    }
    tykohonovParam[i] <- min(max(tykohonovParam[i], 10^-6), 10^8)
  }

  return(tykohonovParam)
}

#' Compute the MLE for a selected region of interest
#'
#' Computes the conditional mle for a region selected based on
#' the selection rule \code{y[selected] > threshold} or
#' \code{y[selected] < -threshold}, and the coordinates which were not
#' selected must violate the selection rule.
#'
#' @param y the observed noraml coordinates
#'
#' @param cov the covariance of \code{y}
#'
#' @param threshold the threshold used in the selection rule.
#' Must be either a scalar or a numeric vector of size \code{length(y)}
#'
#' @param coordinates an optional matrix of the coordinates of the observed
#' vector. This is only relevant if \code{y} corresponds to a spatial
#' or temporal observation
#'
#' @param selected an optional boolean vector, with \code{TRUE} coordinates
#' corresponding to coordinates of \code{y} that were selected.
#'
#' @param projected an optional fixed value that \code{mean(mu)} must equal.
#' Can be used to construct profile likelihood post-selection confidence
#' intervals.
#'
#' @param tykohonovParam an optional penalty value for the tykohonov
#' regularizer. This is an (inferior) altenative to specifying a
#' \code{tykohonovSlack}.
#'
#' @param tykohonovSlack the estimation routine uses first differences
#' Tykohonov regularization to estimate the mean of the selected region.
#' \code{tykohonovSlack} speficies the allowed deviation from the observed
#' first order differences. The description for details
#'
#' @param stepSizeCoef step size coefficients for stochastic gradient step.
#' Best left unchanged.
#'
#' @param stepRate the rate at which the stochastic gradient steps size should
#' decrease as a function of the number of steps already taken.
#'
#' @param sampPerIter the number of slice MH samples to take for computing the
#' stochastic gradient estimate in each stochastic gradient step
#'
#' @param delay the number of iterations to wait before starting to decrease
#' the stochastic gradient step size
#'
#' @param maxiter the number of stochastic gradient steps to take
#'
#' @param assumeConvergence after how many gradient steps should we assume
#' convergence? The final MLE estimate will be the average of the last
#' \code{maxiter - assumeConvergence} estimates
#'
#' @param nsamp the number of samples to take from the estimated
#' post-selection distribution
#'
#' @param init initial value for the mean estimate
#'
#' @param imputeBoundary the boundary imputation method to use. See
#' description for details
#'
#' @import Matrix
#' @import progress
#' @export
roiMLE <- function(y, cov, threshold,
                   coordinates = NULL,
                   selected = NULL,
                   projected = NULL,
                   tykohonovParam = NULL,
                   tykohonovSlack = 1,
                   stepSizeCoef = 0.15,
                   stepRate = 0.65,
                   sampPerIter = 40,
                   delay = 100,
                   maxiter = 2000,
                   assumeConvergence = 1500,
                   nsamp = 0,
                   init = NULL,
                   progress = FALSE,
                   imputeBoundary = c("smooth", "neighbors", "none", "mean")) {
  maxiter <- max(maxiter, assumeConvergence + length(y) + 1)
  # Basic checks and preliminaries ---------
  if(!(length(threshold) %in% c(1, length(y)))) {
    stop("threshold must be either: a scalar, or a numeric vector of size
         length(y).")
  }

  if(length(threshold) == 1) {
    threshold <- rep(threshold, length(y))
  }

  if(is.null(selected)){
    selected <- abs(y) > threshold
  }
  nselected <- sum(selected)
  p <- length(y)
  s <- sum(selected)
  if(!all(sign(y[selected]) == sign(y[selected][1]))) stop("Signs of active set must be identical")

  # Setting up boundary imputation ----------
  if(length(imputeBoundary) > 1) imputeBoundary <- imputeBoundary[1]
  if(imputeBoundary == "neighbors") {
    if(is.null(coordinates)) {
      imputeBoundary <- "mean"
    } else {
      unselected <- which(!selected)
      distances <- as.matrix(dist(coordinates))
      diag(distances) <- Inf
      neighbors <- cbind(unselected, apply(distances[unselected, , drop = FALSE], 1, function(x) which(selected)[which.min(x[selected])[1]]))
    }
  } else if(imputeBoundary == "smooth") {
    if(!is.null(coordinates)) {
      distances <- as.matrix(dist(coordinates), method = "manhattan")
      distances <- distances[!selected, selected]
      SMOOTH_SD <- 1.5
      smooth_weights <- dnorm(distances, sd = SMOOTH_SD)
      if(sum(!selected) == 1) {
        smooth_weights <- smooth_weights / sum(smooth_weights)
      } else {
        smooth_weights <- apply(smooth_weights, 1, function(x) x / sum(x)) %>% t()
      }
    }
  }

  # Setting-up Tykohonov regularization --------------
  if(is.infinite(tykohonovSlack) | sum(selected) <= 1) {
    tykohonovParam <- rep(0, 2)
  } else if(is.null(tykohonovParam)) {
    tykohonovParam <- c(1, 0)
  } else if(any(tykohonovParam < 0)) {
    tykohonovParam <- c(1, 0)
  } else if(is.null(coordinates)) {
    tykohonovSlack <- 2
    tykohonovParam <- c(1, 0)
  } else if(length(tykohonovParam) == 1) {
    tykohonovParam <- rep(tykohonovParam, 2)
  } else {
    tykohonovParam <- tykohonovParam[1:2]
  }

  if(any(tykohonovParam > 0)) {
    tykMat <- computeTykohonov(selected, coordinates)
    firstDiff <- tykMat$firstDiff
    firstDiff <- as.matrix(Matrix::t(firstDiff) %*% firstDiff)
    secondDiff <- tykMat$secondDiff
    secondDiff <- as.matrix(Matrix::t(secondDiff) %*% secondDiff)
  }

  if(!is.null(init)) {
    if(length(init) != length(y)) stop("If init is not NULL, then its length must match the length of y.")
    mu <- init
  } else {
    mu <- y
    mu[!selected] <- 0
  }

  sds <- sqrt(diag(cov))
  vars <- sds^2
  invcov <- solve(cov)
  condSigma <- 1 / diag(invcov)
  suffStat <- as.numeric(invcov %*% y)
  a <- -threshold
  b <- threshold
  signs <- sign(y)
  obsmean <- mean(y[selected])

  # Setting up projection -----------------
  # If a projected gradient method is used then initalization must be from
  # a parameter which satisfies the constraint.
  if(!is.null(projected)) {
    mu <- mu * sqrt(vars)
    mu[selected] <- rep(projected, sum(selected))
    mu <- mu / sqrt(vars)
  }

  estimates <- matrix(nrow = maxiter, ncol = p)
  estimates[1, ] <- mu

  # Initializing Tykohonov Penalization Parameters
  if(any(tykohonovParam > 0)) {
    yselected <- y[selected]
    obsDiff <- rep(NA, 2)
    obsDiff[1] <- as.numeric(crossprod(yselected, crossprod(firstDiff, yselected)))
    obsDiff[2] <- as.numeric(crossprod(yselected, crossprod(secondDiff, yselected)))
  }

  if(!is.null(projected)) {
    slackAdjusted <- FALSE
  } else {
    tykohonovSlack <- tykohonovSlack * mean(mu[selected]) / obsmean
    slackAdjusted <- TRUE
  }
  tykohonovParam <- pmax(tykohonovParam, 10^-3)

  ############################
  # VARIABLES FOR NEW SAMPLER
  currentSamp <- y
  sampsig <- cov
  sampsig[selected, !selected] <- -sampsig[selected, !selected]
  sampsig[!selected, selected] <- -sampsig[!selected, selected]
  chol <- svd(sampsig)
  chol <- chol$u %*% diag(sqrt(chol$d)) %*% t(chol$v)
  signs <- sign(y)
  lth <- a
  uth <- b
  lth[!selected] <- -lth[!selected]
  uth[!selected] <- -uth[!selected]
  sampMat <- matrix(0.0, nrow = sampPerIter, ncol = length(y))
  gradNorm <- diag(invcov)
  #############################

  if(progress) pb <- txtProgressBar(min = 0, max = maxiter, style = 3)
  for(i in 2:maxiter) {
    # SLICE SAMPLING!
    sampInit <- currentSamp
    sampInit[!selected] <- -sampInit[!selected]
    sampmu <- mu
    sampmu[!selected] <- -sampmu[!selected]
    sampInit <- sampInit - sampmu

    sliceSamplerRcpp(sampMat = sampMat, samp = sampInit,
                     chol = chol,
                     lth = lth - sampmu, uth = uth - sampmu)
    sampMat <- t(t(sampMat) + sampmu)
    sampMat[, !selected] <- -sampMat[, !selected]
    samp <- colMeans(sampMat)

    # Computing gradient ---------------
    # if(i == assumeConvergence) {
    #   stepSizeCoef <- 0
    # }
    condExp <- as.numeric(invcov %*% samp)
    if(tykohonovParam[1] > 0 & sum(selected) > 1) {
      firstGrad <- - as.numeric(firstDiff %*% mu[selected]) * tykohonovParam[1]
    } else {
      firstGrad <- 0
    }
    if(tykohonovParam[2] > 0 & sum(selected) > 1) {
      secondGrad <- - as.numeric(secondDiff %*% mu[selected]) * tykohonovParam[2]
    } else {
      secondGrad <- 0
    }
    barrierGrad <- rep(0, length(y))
    barrierGrad[selected] <- firstGrad + secondGrad

    # Computing gradient and projecting if necessary
    # The projection of the gradient is simply setting its mean to zero
    rawGradient <- (suffStat - condExp + barrierGrad)
    gradNorm <- 0.995 * gradNorm + 0.005 * rawGradient^2
    if(i == 2) {
      MOM_CONST <- 0.2
      momentum_grad <- rawGradient * (1 - MOM_CONST) * 2
    } else {
      momentum_grad <- momentum_grad * MOM_CONST + rawGradient * (1 - MOM_CONST)
    }
    effective_step_coef <- stepSizeCoef / max(i - delay, 1)^stepRate
    gradient <-  momentum_grad / sqrt(gradNorm) * effective_step_coef
    #### TEST DIAGNOSTICS -------
    # cat("gradient mean: ", mean(gradient) / effective_step_coef, "\n")
    # cat("gradient momentum: ", momentum_grad, "\n")
    # cat("gradient (sqrt) norm: ", sqrt(gradNorm), "\n")
    ###########################
    gradient[!selected] <- 0
    gradsign <- sign(gradient)
    gradient <- pmin(abs(gradient), 0.1 * sqrt(vars)) * gradsign
    if(!is.null(projected)) {
      # gradient <- gradient * sqrt(vars)
      gradMean <- mean(gradient[selected])
      gradient[selected] <- gradient[selected] - gradMean
      # gradient <- gradient / sqrt(vars)
    }

    # Updating estimate. The error thing is to make sure we didn't accidently
    # cross the barrier. There might be a better way to do this.
    mu <- mu + gradient

    # Boundary imputation
    if(imputeBoundary != "none") {
      if(imputeBoundary == "mean") {
        mu[!selected] <- mean(mu[selected])
      } else if(imputeBoundary == "neighbors") {
        mu[neighbors[, 1]] <- mu[neighbors[, 2]] * abs(y[neighbors[, 1]] / y[neighbors[, 2]])
      } else if(imputeBoundary == "smooth") {
        mu[!selected] <- smooth_weights %*% mu[selected]
      }
    }

    if(is.null(projected)) {
      mu <- pmax(0, mu * sign(y)) * sign(y)
      mu <- pmin(abs(mu), abs(y)) * sign(y)
    }

    # Updating Tykohonov Params
    if(i > assumeConvergence / 3 &
       !slackAdjusted &
       is.null(projected) &
       sum(selected) > 1) {
      tykohonovSlack <- pmax(tykohonovSlack * mean(mu[selected]) / obsmean, 10^-3)
      slackAdjusted <- TRUE
    }

    if(any(tykohonovParam > 0) & sum(selected) > 1) {
      tykohonovParam <- adjustTykohonov(obsDiff, obsmean, mu, selected,
                                        firstDiff, secondDiff,
                                        tykohonovSlack, tykohonovParam)
    }
    estimates[i, ] <- mu
    if(progress) setTxtProgressBar(pb, i)
  }
  if(progress) close(pb)
  #cat("\n")

  # Sampling from estimated mean --------
  sampmu <- colMeans(estimates[assumeConvergence:maxiter, ])
  if(nsamp > 0) {
    # browser()
    initsamp <- y
    initsamp[!selected] <- -initsamp[!selected]
    sampmu[!selected] <- -sampmu[!selected]
    outSamples <- sliceRcppInner(nsamp, initsamp, sampmu, chol, lth, uth)
    outSamples[, !selected] <- -outSamples[, !selected]
  } else {
    outSamples <- NULL
  }

  # Unnormalizing estimates and samples --------------------------
  # for(i in 1:length(y)) {
  #   if(nsamp > 0) {
  #     outSamples[, i] <- outSamples[, i] * sqrt(vars[i])
  #   }
  #   estimates[, i] <- estimates[, i] * sqrt(vars[i])
  # }

  conditional <- colMeans(estimates[assumeConvergence:maxiter, ])
  return(list(sample = outSamples,
              estimates = estimates,
              conditional = conditional))
}
