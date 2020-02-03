library(dplyr)
library(reshape2)
library(ggplot2)
library(selectiveROI)

#Generating Data ---------------
set.seed(9)
rho <- 0.9
p <- 100
sigma <- rho^as.matrix(dist(1:p))
mu <- rep(0, p) #
mu <- predict(smooth.spline(rnorm(p, sd = 4), df = p/5))$y
y <- as.numeric(mvtnorm::rmvnorm(1, mean = mu, sigma = sigma))
threshold <- 1.5 / sqrt(diag(sigma))
select <- abs(y) > 1.5

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

# sampling -----
# probability under true param
# Computing MLE using old MLE function -----
system.time(newmle1 <- roiMLE(y, sigma, threshold, compute = "mle",
                 coordinates = matrix(1:length(y), ncol = 1),
                 selected = select,
                 regularization_slack = 1))
system.time(newmle2 <- roiMLE(y, sigma, threshold, compute = "mle",
                             coordinates = matrix(1:length(y), ncol = 1),
                             selected = select,
                             regularization_slack = 2))

system.time(lcifit <- roiMLE(y, sigma, threshold, compute = "lower-CI",
                             coordinates = matrix(1:length(y), ncol = 1),
                             selected = select,
                             regularization_slack = 2))
lci = lcifit$ci_lim_path %>% tail(1)
system.time(newmle <- roiMLE(y, sigma, threshold, compute = "upper-CI",
                             coordinates = matrix(1:length(y), ncol = 1),
                             selected = select,
                             regularization_slack = 1))
uci = newmle$ci_lim_path %>% tail(1)

ylim_l = min(c(lci, -threshold)) - 1
ylim_u = max(c(uci, threshold )) + 1
plot(keep, y, type = "l", col = "red", ylim = c(ylim_l, ylim_u))
lines(keep, mu, col = "black")
abline(h = c(-threshold, threshold), col = "grey", lty = 2)
abline(v = c(min(keep - 0.1), max(keep + 0.1)), lty = 2)
lines(keep, newmle1$conditional, col = "dark green")
lines(keep, newmle2$conditional, col = "orange")
lines(keep, rep(lci, length(keep)), col = "blue")
lines(keep, rep(uci, length(keep)), col = "blue")
lines(keep, rep(mean(mu[select]), length(keep)), col = "black", lty = 2)
lines(keep, rep(mean(y[select]), length(keep)), col = "red", lty = 2)
abline(h = 0, lty = 3)




