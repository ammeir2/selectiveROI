# verifies that the ROI selection rule is valid.
# returns the seleced coordinates or throws an error
#
# y - the observed vector
# threshold - a vector of size length(y)
# selected - an optional boolean vector of size length(y)
verifySelection <- function(y, threshold, selected = NULL) {
  # Identifying selected coordinates
  anyPos <- any(y >= threshold)
  anyNeg <- any(y <= -threshold)
  if(is.null(selected)) {
    if(anyPos & anyNeg) {
      stop("Cant identify selected coordinates, both positive crossing
           and negative crossing observations in observed vector!")
    } else if(anyPos) {
      selected <- y >= threshold
    } else if(anyNeg) {
      selected <- y <= -threshold
    } else {
      stop("No selected coordinates!")
    }
    return(selected)
  }

  # Verifying selection event
  if(!any(selected)) {
    stop("No selected coordinates?!")
  }

  if(anyPos) {
    if(!all(y[selected] >= threshold[selected])) {
      stop("not all selected coordinates cross the threshold")
    }
    if(!all(y[!selected] < threshold[!selected])) {
      stop("not all of the excluded coordinates are below threshold!")
    }
  }

  if(anyNeg) {
    if(!all(y[selected] <= -threshold[selected])) {
      stop("not all selected coordinates cross the threshold")
    }
    if(!all(y[!selected] > -threshold[!selected])) {
      stop("not all of the excluded coordinates are above threshold!")
    }
  }

  return(selected)
}

#' A slice sampler for sampling from the ROI selection event
#'
#' Selection event is assumed to be, {selected > threshold} or
#' {selected < -threshold} (with coordinates which were not selected not
#' satisfying the selection event)
#'
#' @param nsamp number of samples to be taken
#'
#' @param y the observed normal vector, or an initial sample point that
#' satisfies the selection event
#'
#' @param mu the mean vector of the (untruncated) normal vector
#'
#' @param sigma the covariance of the (untruncated) normal vector
#'
#' @param threshold a vector of size length(y), the selection threshold
#'
#' @param selected which coordinates were selected? if NULL function
#' will attempt to figure it out on its own but my fail if there are
#' both observations that are larger than \code{threshold} and observations
#' that are smaller than \code{-threshold}.
#'
#' @importFrom Rcpp evalCpp
#' @import magrittr
#' @useDynLib selectiveROI
#' @export
sliceSamplerROI <- function(nsamp, y, mu, sigma, threshold, selected = NULL) {
  # verifying threshold
  if(length(threshold) == 1) {
    threshold <- rep(threshold, length(y))
  } else if(length(threshold) != length(y)) {
    stop("threshold must either be a real number or a vector of size length(y)!")
  }

  # Veriftiny selection
  selected <- verifySelection(y, threshold, selected)

  # Transforming to same side selection event
  signs <- sign(y)
  y[selected] <- signs[selected] * y[selected]
  y[!selected] <- -signs[!selected] * y[!selected]
  mu[selected] <- signs[selected] * mu[selected]
  mu[!selected] <- -signs[!selected] * mu[!selected]
  threshold[!selected] <- -threshold[!selected]
  sigma[!selected, selected] <- -sigma[!selected, selected]
  sigma[selected, !selected] <- -sigma[selected, !selected]
  chol <- svd(sigma)
  chol <- chol$u %*% diag(sqrt(chol$d)) %*% t(chol$v)
  samples <- sliceRcppInner(nsamp, init = y, mu = mu,
                            chol = chol,
                            lth = -threshold, uth = threshold)
  for(j in 1:ncol(samples)) {
    if(selected[j]) {
      samples[, j] <- samples[, j] * signs[j]
    } else {
      samples[, j] <- samples[, j] * signs[j] * -1
    }
  }

  return(samples)
}

# The inner part of the slice sampling routine
#
# nsamp - the number of samples
# init - the starting point of the sampler
# mu - the mean of the normal vector
# chol - square root of the covariance
# lth - the lower selection threshold
# uth - the upper selection threshold
sliceRcppInner <- function(nsamp, init, mu, chol, lth, uth) {
  init <- init - mu
  samp <- init
  sampMat <- matrix(0.0, nrow = nsamp, ncol = length(init))
  sliceSamplerRcpp(sampMat = sampMat, samp = samp,
                   chol = chol, lth = lth - mu, uth = uth - mu)
  sampMat <- apply(sampMat, 1, function(x) x + mu) %>% t()
  return(sampMat)
}

#' A function for computing the probability of positive/negative selections
#'
#' This function uses the mvtnorm package to compute the probability of
#' positive and negative selection events of selective ROI problems.
#' The selection event is assumed to be, {selected > threshold} or
#' {selected < -threshold} (with coordinates which were not selected not
#' satisfying the selection event)
#'
#' @param y the observed normal vector, or an initial sample point that
#' satisfies the selection event
#'
#' @param mu the mean vector of the (untruncated) normal vector
#'
#' @param sigma the covariance of the (untruncated) normal vector
#'
#' @param threshold a vector of size length(y), the selection threshold
#'
#' @param selected which coordinates were selected? if NULL function
#' will attempt to figure it out on its own but my fail if there are
#' both observations that are larger than \code{threshold} and observations
#' that are smaller than \code{-threshold}.
#'
#' @importFrom mvtnorm pmvnorm
#' @export
mvtnormSelectionProbs <- function(y, mu, sigma, threshold, selected = NULL) {
  # verifying threshold
  if(length(threshold) == 1) {
    threshold <- rep(threshold, length(y))
  } else if(length(threshold) != length(y)) {
    stop("threshold must either be a real number or
         a vector of size length(y)!")
  }

  # Verifying selection
  selected <- verifySelection(y, threshold, selected)

  # Transforming to same side selection event
  signs <- sign(y)
  mu[selected] <- signs[selected] * mu[selected]
  mu[!selected] <- -signs[!selected] * mu[!selected]
  sigma[!selected, selected] <- -sigma[!selected, selected]
  sigma[selected, !selected] <- -sigma[selected, !selected]
  threshold[!selected] <- threshold[!selected]

  # Computing probabilities
  posProb <- pmvnorm(lower = threshold, upper = Inf,
                     mean = mu, sigma = sigma)
  negProb <- pmvnorm(lower = -Inf, upper = -threshold,
                     mean = mu, sigma = sigma)
  posProb <- posProb / (posProb + negProb)
  negProb <- 1 - posProb

  # Transforming back to original signs
  if(sign(y[selected][1]) == 1) {
    result <- c(positive = posProb, negProb = negProb)
  } else {
    result <- c(positive = negProb, negProb = posProb)
  }
  return(result)
}
