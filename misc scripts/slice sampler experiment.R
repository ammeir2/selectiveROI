library(dplyr)
library(reshape2)
library(ggplot2)
library(selectiveROI)

#Generating Data ---------------
set.seed(2)
rho <- 0.9
p <- 100
sigma <- rho^as.matrix(dist(1:p))
mu <- rep(0, p) #
mu <- predict(smooth.spline(rnorm(p, sd = 2), df = p/5))$y
y <- as.numeric(mvtnorm::rmvnorm(1, mean = mu, sigma = sigma))
plot(y, type = "l", col = "red", ylim = c(min(c(y, mu)), max(c(y, mu))))
lines(mu, col = "black")
threshold <- 1.5 / sqrt(diag(sigma))
select <- abs(y) > 1.5
abline(h = c(-threshold, threshold), col = "grey", lty = 2)

# Identifying clusters --------
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
abline(v = c(min(keep - 0.1), max(keep + 0.1)), lty = 2)

# sampling -----
# probability under true param
length(y)
system.time(sample <- sliceSamplerROI(5000, y, mu, sigma, threshold, selected = select))
mvtProbs <- mvtnormSelectionProbs(y, mu, sigma, threshold, selected = select)
mean(sample[, 2] > 0)
mvtProbs

# probability under the null
zeros <- rep(0, length(y))
nullSample <- sliceSamplerROI(5000, y, zeros, sigma, threshold, selected = select)
mean(nullSample[, 2] > 0)

# Computing MLE using old MLE function -----
system.time(newmle <- roiMLE(y, sigma, threshold,
                 coordinates = matrix(1:length(y), ncol = 1),
                 selected = select,
                 stepSizeCoef = 0.15, stepRate = 0.65,
                 delay = 100, maxiter = 4010,
                 trimSample = 50,
                 assumeConvergence = 4000,
                 probMethod = "all",
                 imputeBoundary = "neighbors"))
lines(keep, newmle$conditional, col = "dark green")

system.time(oldmle <- optimizeSelected(y, sigma, threshold,
                           coordinates = matrix(1:length(y), ncol = 1),
                           selected = select,
                           stepSizeCoef = 0.25, stepRate = 0.65,
                           delay = 100, maxiter = 2000,
                           assumeConvergence = 1950,
                           probMethod = "all",
                           imputeBoundary = "neighbors"))
lines(keep, oldmle$conditional, col = "blue")
mean((mu[select] - oldmle$conditional[select])^2)
mean((mu[select] - newmle$conditional[select])^2)




