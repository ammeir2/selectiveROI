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
#' @param mean_weights the weights to use for the contrast to be computed, if
#' not specified then equal weights will be given to all selected coordinates.
#' \code{mean_weights} do not have to sum to one. If specified, then the length
#' of the vector must either equal \code{sum(selected)} (the number of selected
#' coordinates), or \code{length(y)}. Either way, all coordinates which were not
#' selected must be given a weight of zero.
#'
#' @param projected an optional fixed value that \code{mean(mu)} must equal.
#' Can be used to construct profile likelihood post-selection confidence
#' intervals.
#'
#' @param regularization_param an optional penalty value for the tykohonov
#' regularizer. This is an (inferior) altenative to specifying a
#' \code{regularization_slack}.
#'
#' @param regularization_slack the estimation routine uses first differences
#' Tykohonov regularization to estimate the mean of the selected region.
#' \code{regularization_slack} speficies the allowed deviation from the observed
#' first order differences. The description for details
#'
#' @param init initial value for the mean estimate
#'
#' @param progress whether to display a bar describing the progress
#' of the gradient algorithm.
#'
#' @param sampling_control a list with control parameters for sampling
#' from the estimated distribution.
#'
#' @param mle_contorl a list of parameters to be used when computed the
#' conditional MLE.
#'
#' @import Matrix
#' @import progress
#' @export
roiMLE <- function(y, cov, threshold,
                   compute = c("mle", "lower-CI", "upper-CI"),
                   ci_alpha = 0.025,
                   coordinates = NULL,
                   selected = NULL,
                   mean_weights = NULL,
                   projected = NULL,
                   regularization_param = NULL,
                   regularization_slack = 1,
                   init = NULL,
                   progress = FALSE,
                   sampling_control = roi_sampling_control(),
                   mle_control = roi_mle_control()) {
  list2env(mle_control, environment())
  grad_iterations <- max(grad_iterations, assume_convergence + length(y) + 1)
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
  if(length(impute_boundary) > 1) impute_boundary <- impute_boundary[1]
  if(impute_boundary == "neighbors") {
    if(is.null(coordinates)) {
      stop("Must specify coordinates for impute_boundary == 'neighbors'")
    } else {
      unselected <- which(!selected)
      distances <- as.matrix(dist(coordinates))
      diag(distances) <- Inf
      neighbors <- cbind(unselected, apply(distances[unselected, , drop = FALSE], 1, function(x) which(selected)[which.min(x[selected])[1]]))
    }
  } else if(impute_boundary == "smooth") {
    if(!is.null(coordinates)) {
      distances <- as.matrix(dist(coordinates), method = "manhattan")
      distances <- distances[!selected, selected]
      SMOOTH_SD <- 0.25
      smooth_weights <- dnorm(distances, sd = SMOOTH_SD)
      if(sum(!selected) == 1) {
        smooth_weights <- smooth_weights / sum(smooth_weights)
      } else {
        smooth_weights <- apply(smooth_weights, 1, function(x) x / sum(x)) %>% t()
      }
    }
  }

  # Setting-up Tykohonov regularization --------------
  # If only one coordinate is selected then we don't need to regularize
  if(is.infinite(regularization_slack) | sum(selected) <= 1) {
    regularization_param <- rep(0, 2)
  } else if(is.null(regularization_param)) {
    regularization_param <- c(1, 0)
  } else if(any(regularization_param < 0)) {
    stop("Regularization parameters must be positive")
  } else if(is.null(coordinates)) {
    regularization_slack <- 2
    regularization_param <- c(1, 0)
  } else if(length(regularization_param) == 1) {
    regularization_param <- rep(regularization_param, 2)
  } else {
    regularization_param <- regularization_param[1:2]
  }

  if(any(regularization_param > 0)) {
    tykMat <- computeTykohonov(selected, coordinates)
    firstDiff <- tykMat$firstDiff
    firstDiff <- as.matrix(Matrix::t(firstDiff) %*% firstDiff)
    secondDiff <- tykMat$secondDiff
    secondDiff <- as.matrix(Matrix::t(secondDiff) %*% secondDiff)
  }

  if(!is.null(init)) {
    if(length(init) != length(y)) {
      stop("If init is not NULL, then its length must match the length of y.")
    }
    mu <- init
  } else {
    mu <- y
    mu[!selected] <- 0
  }

  sds <- sqrt(diag(cov))
  vars <- diag(cov)
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
    mu[selected] <- rep(projected, sum(selected))
  }

  estimates <- matrix(nrow = grad_iterations, ncol = p)
  estimates[1, ] <- mu

  # Initializing Tykohonov Penalization Parameters
  if(any(regularization_param > 0)) {
    yselected <- y[selected]
    obsDiff <- rep(NA, 2)
    obsDiff[1] <- as.numeric(crossprod(yselected, crossprod(firstDiff, yselected)))
    obsDiff[2] <- as.numeric(crossprod(yselected, crossprod(secondDiff, yselected)))
  }

  if(!is.null(projected)) {
    slackAdjusted <- FALSE
    effective_slack <- regularization_slack
  } else {
    effective_slack <- regularization_slack * mean(mu[selected]) / obsmean
    slackAdjusted <- TRUE
  }
  regularization_param <- pmax(regularization_param, 10^-3)

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
  sampMat <- matrix(0.0, nrow = samp_per_iter, ncol = length(y))
  gradNorm <- diag(invcov)
  #############################

  ### New projection method ########################
  weight_stop_message <-"length(mean_weights) must either be equal to length(y) or sum(selected),
  coordinates which were not selected are must have a weight of zero."
  if(is.null(mean_weights)) {
    mean_weights <- rep(1, sum(selected)) / sum(selected)
  } else {
    if(length(mean_weights) == length(y)) {
      if(any(mean_weights[!selected] != 0)) {
        stop(weight_stop_message)
      }
      mean_weights <- mean_weights[selected]
    } else if(length(mean_weights) != sum(selected)) {
      stop(weight_stop_message)
    }
  }

  compute <- compute[1]
  if(!is.null(projected) | compute != "mle") {
    sub_invcov <- invcov[selected, selected]
    sub_invcov <- sub_invcov*0.5  + 0.5*diag(diag(sub_invcov))
    mahal_vec <- as.numeric(sub_invcov %*% mean_weights)
    mahal_const <- sum(mahal_vec * mean_weights)
  }
  ##########################################

  # Saving some samples
  samples_for_init <- matrix(NA, nrow = ceiling((grad_iterations - assume_convergence) / 10), ncol = length(y))
  init_mat_row <- 1

  # What are we computing? ---------
  if(compute %in% c("lower-CI", "upper-CI")) {
    mean_sd <- as.numeric(t(mean_weights) %*% cov[selected, selected] %*% mean_weights)
    obs_mean <- sum(y[selected] * mean_weights)
    RB_const <- 2 / dnorm(2 * ci_alpha, sd = mean_sd)
  }

  if(compute == "lower-CI") { # Initializing from a conservative limit
    if(obs_mean > 0) {
      ci_lim_level = ci_alpha^2
    } else {
      ci_lim_level = ci_alpha
    }
    ci_lim <- obs_mean - qnorm(1 - ci_lim_level, sd = mean_sd)
    effective_slack <- regularization_slack * abs(ci_lim) / abs(obs_mean)
    mu <- project_vector_mahalanobis(y, ci_lim, selected, mean_weights,
                                     mahal_const, mahal_vec)
    ci_lim_path <- numeric(grad_iterations - 1)
    samp_means <- numeric(grad_iterations - 2)
    ci_lim_path[1] <- ci_lim
  } else if(compute == "upper-CI") {
    if(obs_mean < 0) {
      ci_lim_level = ci_alpha^2
    } else {
      ci_lim_level = ci_alpha
    }
    ci_lim <- obs_mean + qnorm(1 - ci_alpha^2, sd = mean_sd)
    effective_slack <- regularization_slack * abs(ci_lim) / abs(obs_mean)
    mu <- project_vector_mahalanobis(y, ci_lim, selected, mean_weights,
                                     mahal_const, mahal_vec)
    ci_lim_path <- numeric(grad_iterations - 1)
    samp_means <- numeric(grad_iterations - 2)
    ci_lim_path[1] <- ci_lim
  } else if(compute != "mle") {
    stop("compute must be either, lower-CI, upper-CI, or mle")
  } else {
    ci_lim_path <- NULL
  }

  if(progress) pb <- txtProgressBar(min = 0, max = grad_iterations, style = 3)
  for(i in 2:grad_iterations) {
    # Projecting to new CI limit
    if(i > 2 & compute != "mle") {
      STEP_LIM <- 0.05
      if(compute == "lower-CI") {
        mult_const <- RB_const * (obs_mean - ci_lim)
        if(samp_mean < obs_mean) {
          ci_lim <- ci_lim + min(mult_const * ci_alpha / sqrt(i - 2), STEP_LIM * mean_sd)
        } else {
          ci_lim <- ci_lim - min(mult_const * (1 - ci_alpha) / sqrt(i - 2), STEP_LIM * mean_sd)
        }
        effective_slack <- regularization_slack * min(abs(ci_lim) / abs(obs_mean), abs(obs_mean - mean_sd * qnorm(1 - ci_alpha)) / abs(obs_mean))
        #print(effective_slack)
      } else if(compute == "upper-CI") {
        mult_const <- RB_const * (ci_lim - obs_mean)
        if(samp_mean > obs_mean) {
          ci_lim <- ci_lim - min(mult_const * ci_alpha / (i - 2), STEP_LIM * mean_sd)
        } else {
          ci_lim <- ci_lim + min(mult_const * (1 - ci_alpha) / (i - 2), STEP_LIM * mean_sd)
        }
        effective_slack <- regularization_slack * min(abs(ci_lim) / abs(obs_mean), abs(obs_mean + mean_sd * qnorm(1 - ci_alpha)) / abs(obs_mean))
      }
      ci_lim_path[i - 1] <- ci_lim
      samp_means[i - 2] <- samp_mean
      mu <- project_vector_mahalanobis(mu, ci_lim, selected, mean_weights,
                                       mahal_const, mahal_vec)
    }

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
    if(compute != "mle") {
      samp_mean <- sum(sampMat[nrow(sampMat), selected] * mean_weights)
    }

    if(i > assume_convergence & i %% 10 == 0) {
      # print(c(nrow(samples_for_init), init_mat_row))
       samples_for_init[init_mat_row, ] <- sampMat[nrow(sampMat), ]
      init_mat_row <- init_mat_row + 1
    }

    condExp <- as.numeric(invcov %*% samp)
    if(regularization_param[1] > 0 & sum(selected) > 1) {
      firstGrad <- - as.numeric(firstDiff %*% mu[selected]) * regularization_param[1]
    } else {
      firstGrad <- 0
    }
    if(regularization_param[2] > 0 & sum(selected) > 1) {
      secondGrad <- - as.numeric(secondDiff %*% mu[selected]) * regularization_param[2]
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
    effective_step_coef <- step_size_coef / max(i - grad_delay, 1)^step_rate
    gradient <-  momentum_grad / sqrt(gradNorm) * effective_step_coef
    #### TEST DIAGNOSTICS -------
    # cat("gradient mean: ", mean(gradient) / effective_step_coef, "\n")
    # cat("gradient momentum: ", momentum_grad, "\n")
    # cat("gradient (sqrt) norm: ", sqrt(gradNorm), "\n")
    ###########################
    gradient[!selected] <- 0
    gradsign <- sign(gradient)
    gradient <- pmin(abs(gradient), 0.1 * sqrt(vars)) * gradsign
    if(!is.null(projected) | compute != "mle") {
      gradient <- project_vector_mahalanobis(gradient, 0, selected, mean_weights,
                                             mahal_const, mahal_vec)
    }

    # Updating estimate. The error thing is to make sure we didn't accidently
    # cross the barrier. There might be a better way to do this.
    mu <- mu + gradient

    # Boundary imputation
    if(impute_boundary != "none") {
      if(impute_boundary == "mean") {
        mu[!selected] <- mean(mu[selected])
      } else if(impute_boundary == "neighbors") {
        mu[neighbors[, 1]] <- mu[neighbors[, 2]] * abs(y[neighbors[, 1]] / y[neighbors[, 2]])
      } else if(impute_boundary == "smooth") {
        if(i == 2) {
          neighbor_ratio <- as.numeric(abs(y[!selected]) / abs(smooth_weights %*% y[selected]))
        }
        mu[!selected] <- smooth_weights %*% mu[selected] * neighbor_ratio
      }
    }

    if(compute == "mle" & is.null(projected)) {
      mu <- pmax(0, mu * sign(y)) * sign(y)
      mu <- pmin(abs(mu), abs(y)) * sign(y)
    }

    # Updating Tykohonov Params
    if(i > assume_convergence / 3 &
       !slackAdjusted &
       is.null(projected) &
       sum(selected) > 1 &
       compute == "mle") {
      effective_slack <- pmax(regularization_slack * mean(mu[selected]) / obsmean, 10^-3)
      slackAdjusted <- TRUE
    }

    if(any(regularization_param > 0) & sum(selected) > 1) {
      regularization_param <- adjustTykohonov(obsDiff, obsmean, mu, selected,
                                        firstDiff, secondDiff,
                                        effective_slack, regularization_param)
    }
    estimates[i, ] <- mu
    if(progress) setTxtProgressBar(pb, i)
  }
  if(progress) close(pb)
  samples_for_init <- samples_for_init[1:(init_mat_row - 1), ]
  #cat("\n")

  # Sampling from estimated mean --------
  sampmu <- colMeans(estimates[assume_convergence:grad_iterations, ])
  nsamp <- sampling_control[["samp_per_chain"]]
  n_chains <- sampling_control[["n_chains"]]
  burnin <- sampling_control[["burnin"]]
  if(nsamp > 0) {
    sample_list <- vector(n_chains, mode = "list")
    if(n_chains > 1 & nrow(samples_for_init) > 1) {
      init_replace <- n_chains > nrow(samples_for_init) + 1
      # browser()
      samples_for_init <- samples_for_init[sample.int(nrow(samples_for_init),
                                                       size = n_chains - 1 ,
                                                       replace = init_replace), ]
    }

    for(i in 1:n_chains) {
      # browser()
      if(i == 1 | nrow(samples_for_init) < 2) {
        initsamp <- y
      } else {
        # browser()
        initsamp <- samples_for_init[i - 1, ]
      }
      initsamp[!selected] <- -initsamp[!selected]
      sampmu[!selected] <- -sampmu[!selected]
      sample_list[[i]] <- sliceRcppInner(nsamp, initsamp, sampmu, chol, lth, uth)
      sample_list[[i]][, !selected] <- -sample_list[[i]][, !selected]
      sample_list[[i]] <- sample_list[[i]][(burnin + 1):nrow(sample_list[[i]]), ]
    }
    outSamples <- do.call("rbind", sample_list)
  } else {
    outSamples <- NULL
  }

  conditional <- colMeans(estimates[assume_convergence:grad_iterations, ])
  return(list(sample = outSamples,
              estimates = estimates,
              conditional = conditional,
              ci_lim_path = ci_lim_path))
}

#' Generate a list of parameters for sampling from the estimated conditional distribution
#'
#' A function which generates a list of parameters for use in the roiMLE function,
#' for sampling from the estimated distribution. This can be used for constructing
#' confidence intervals or conducting hypothesis testing.
#'
#' @param samp_per_chain the number of samples to take from the estimated
#' post-selection distribution in each MCMC chain. If 0, then no sampling
#' will be conducted (default).
#'
#' @param burnin the burnin periods for each MCMC chain.
#'
#' @param n_chains the number of MCMC chains to run. The MCMC chains use
#' different intializations. One chain will usually suffice, but it is
#' recommended to run several if there is strong depedence in the data
#' that may slow down the mixing of the MCMC chain.
#'
#' @export
roi_sampling_control <- function(samp_per_chain = 0,
                             burnin = 1000,
                             n_chains = 5) {
  control <- list(samp_per_chain = samp_per_chain,
                  burnin = burnin,
                  n_chains = n_chains)
  return(control)
}

#' Generate a list of parameters controllong the stochastic gradient process of the \code{roiMLE} function
#'
#' @param grad_iterations the number of stochastic gradient steps to take
#'
#' @param step_size_coef step size coefficients for stochastic gradient step.
#' Best left unchanged.
#'
#' @param step_rate the rate at which the stochastic gradient steps size should
#' decrease as a function of the number of steps already taken.
#'
#' @param samp_per_iter the number of slice MH samples to take for computing the
#' stochastic gradient estimate in each stochastic gradient step
#'
#' @param grad_delay the number of iterations to wait before starting to decrease
#' the stochastic gradient step size
#'
#' @param assume_convergence after how many gradient steps should we assume
#' convergence? The final MLE estimate will be the average of the last
#' \code{grad_iterations - assume_convergence} estimates
#'
#' @param impute_boundary the boundary imputation method to use. See
#' description for details
#'
#' @param RB_mult adjusts the Robins-Monroe step sizes when computing
#' profile-likelihood confidence intervals.
#'
#' @import magrittr
#'
#' @export
roi_mle_control <- function(grad_iterations = 2100,
                        step_size_coef = 0.5,
                        step_rate = 0.55,
                        samp_per_iter = 20,
                        grad_delay = NULL,
                        assume_convergence = NULL,
                        impute_boundary = c("smooth", "neighbors", "none", "mean"),
                        RB_mult = 1) {
  if(is.null(grad_delay)) {
    grad_delay <- min(ceiling(grad_iterations / 6), 100)
  }
  if(is.null(assume_convergence)) {
    assume_convergence <- floor(grad_iterations / 3)
  }
  control <- list(grad_iterations = grad_iterations,
                  step_size_coef = step_size_coef,
                  step_rate = step_rate,
                  samp_per_iter = samp_per_iter,
                  grad_delay = grad_delay,
                  assume_convergence = assume_convergence,
                  impute_boundary = impute_boundary[1],
                  RB_mult = RB_mult)
  return(control)
}

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

#' Adjusts tykohonov parameter
adjustTykohonov <- function(obsDiff, obsmean, mu, selected,
                            firstDiff, secondDiff,
                            tykohonvSlack, regularization_param) {
  mu <- mu[selected]
  meanmu <- mean(mu)
  for(i in 1:2) {
    if(regularization_param[i] == 0) next

    if(i == 1) {
      muDiff <- crossprod(mu, crossprod(firstDiff, mu))
    } else {
      muDiff <- crossprod(mu, crossprod(secondDiff, mu))
    }

    # if(i == 1 & any(sign(mu[1]) != sign(mu)) & mean(mu) > 0.5) {
    #   regularization_param[i] <- regularization_param[i] * 1.1
    #   next
    # }

    ratio <- muDiff / (obsDiff[i] * (tykohonvSlack))
    #if(i == 2) print(c(ratio, muDiff, obsDiff[i]))
    if(is.na(ratio)) {
      regularization_param[i] <- Inf
    } else if(ratio > 2) {
      regularization_param[i] <- regularization_param[i] * 1.3
    } else if(ratio > 1) {
      regularization_param[i] <- regularization_param[i] * 1.05
    } else if(ratio < 0.5) {
      regularization_param[i] <- regularization_param[i] * 0.7
    } else if(ratio < 1) {
      regularization_param[i] * 0.95
    }
    regularization_param[i] <- min(max(regularization_param[i], 10^-6), 10^8)
  }

  return(regularization_param)
}

#' A function for projecting a vector such that it will satisfy a linear constraint
#' based on a mahalanobis metric
#'
#' @param x the vector to be projected
#' @param target what the linear function of x should equal after the projection
#' @param selected which coordinates of x were selected
#' @param mean_weights we require that \code{sum(x * mean_weights)} will equal \code{target}
#' @param mahal_const Let \code{A} be a matrix such that our mahalanobis metric is \code{t(x) %*% A %*% x}.
#' Then \code{mahal_const} is \code{t(mean_weights) %*% A %*% mean_weights}
#' @param mahal_vec \code{A %*% mean_weights}
project_vector_mahalanobis <- function(x, target, selected, mean_weights,
                                       mahal_const, mahal_vec) {
  selected_x <- x[selected]
  proj_adjust <- (target - sum(mean_weights * selected_x)) / mahal_const * mahal_vec
  x[selected] <- selected_x + proj_adjust
  return(x)
}

