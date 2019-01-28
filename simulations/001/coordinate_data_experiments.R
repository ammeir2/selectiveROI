library(foreach)
library(magrittr)
library(doParallel)

# simulation parameters ---------
run_sim <- function(config) {
  set.seed(config[["seed"]])

  data <- generate_smooth_coordinate_data(n = 1,
                                          dims = config[["dims"]][[1]],
                                          rho = 0.95,
                                          noise_sd = config[["noise_sd"]],
                                          mu_sd = config[["mu_sd"]],
                                          bandwidth_sd = config[["bandwidth_sd"]],
                                          dist_method = "manhattan")
  threshold <- config[["threshold"]]
  lims <- max(abs(c(data$y[1, ], data$mu[1, ])))
  lims <- c(-lims, lims)
  # plot(data$y[1, ], type = "l", ylim = lims)
  # lines(data$mu[1, ], type = "l", col = "red")
  # abline(h = c(threshold, - threshold), col = "grey", lty = 2)

  results_names <- c("n_selected", "size", "obs_mean", "true_mean", "cond_mean",
                     "l_quant", "u_quant", "covered")
  repetitions <- nrow(data$y)
  result_list <- vector(repetitions, mode = "list")
  for(m in 1:repetitions) {
    y <- data$y[m, ]
    n_coordinates <- length(y)
    mu <- data$mu[m, ]
    coordinates <- data$coordinates
    covariance <- data$covariance
    distances <- data$distances
    pos_selected <- which(y > threshold)
    pos_clusters <- find_clusters(distances, pos_selected,
                                  dist_threshold = config[["distance_threshold"]],
                                  dist_method = "manhattan")
    neg_selected <- which(y < -threshold)
    neg_clusters <- find_clusters(distances, neg_selected,
                                  dist_threshold = config[["distance_threshold"]],
                                  dist_method = "manhattan")
    covariance <- data$covariance
    clusters <- c(neg_clusters, pos_clusters)
    clusters <- clusters[sapply(clusters, length) >= config[["min_size"]]]
    if(length(clusters) == 0) next
    iter_results <- vector(length(clusters), mode = "list")
    for(i in 1:length(clusters)) {
      this_cluster <- clusters[[i]]
      print(sprintf("Processing cluster %s out of %s size %s", i, length(clusters), length(this_cluster)))
      boundary <- (distances[setdiff(1:n_coordinates, this_cluster), this_cluster] <= config[["distance_threshold"]]) %>%
        apply(1, any) %>% which()
      boundary <- setdiff(1:n_coordinates, this_cluster)[boundary]
      selected <- c(rep(TRUE, length(this_cluster)), rep(FALSE, length(boundary)))
      this_cluster <- c(this_cluster, boundary)
      this_covariance <- covariance[this_cluster, this_cluster]
      mle <- roiMLE(y = y[this_cluster],
                    cov = this_covariance, threshold,
                    coordinates = coordinates[this_cluster, , drop = FALSE],
                    selected = selected,
                    projected = NULL,
                    tykohonovParam = NULL,
                    tykohonovSlack = 1,
                    stepSizeCoef = 0.25,
                    stepRate = 0.85,
                    sampPerIter = 20,
                    delay = 250,
                    maxiter = 2000,
                    assumeConvergence = 750,
                    nsamp = 0,
                    imputeBoundary = c("neighbors"),
                    progress = TRUE)
      # if(FALSE) {
      #   plot_coord <- 1
      #   plot(mle$estimates[, plot_coord])
      #   abline(h = mle$conditional[plot_coord], col = "red")
      #   abline(h = y[this_cluster][plot_coord], col = "green")
      #   abline(h = mu[this_cluster][plot_coord], col = "blue")
      #   rbind(y[this_cluster], mle$conditional, mu[this_cluster])
      #   mean((y[this_cluster][selected] - mu[this_cluster][selected])^2)
      #   mean((mle$conditional[selected] - mu[this_cluster][selected])^2)
      # }

      true_mean <- mean(mu[this_cluster][selected])
      obs_mean <- mean(y[this_cluster][selected])
      profile_mle <- roiMLE(y = y[this_cluster],
                            cov = this_covariance, threshold,
                            coordinates = coordinates[this_cluster, , drop = FALSE],
                            selected = selected,
                            projected = true_mean,
                            tykohonovParam = NULL,
                            tykohonovSlack = abs(true_mean / obs_mean)^config[["tyk_exp"]],
                            stepSizeCoef = 0.25,
                            stepRate = 0.85,
                            sampPerIter = 20,
                            delay = 100,
                            maxiter = 1000,
                            assumeConvergence = 200,
                            nsamp = 10^4 * 2,
                            imputeBoundary = c("neighbors"),
                            progress = TRUE)
      samp <- rowMeans(profile_mle$sample[, selected])
      quants <- quantile(samp, probs = c(0.025, 0.975))
      covered <- quants[1] < true_mean & quants[2] > true_mean
      iter_results[[i]] <- c(n_selected = sum(selected),
                             size = length(this_cluster),
                             true_mean = true_mean,
                             obs_mean = obs_mean,
                             mle = mean(mle$conditional[selected]),
                             q025 = quants[1],
                             q975 = quants[2],
                             covered = covered)
      print(iter_results[[i]])
    }
    result_list[[m]] <- do.call("rbind", iter_results)
    lapply(result_list[1:m], colMeans) %>% do.call("rbind", .) %>% colMeans() %>% print()
  }
  return(list(parameters = config,
              results = result_list))
}

configurations <- expand.grid(seed = 1:100,
                              dims = list(c(200, 1), c(14, 14)),
                              threshold = c(1.65, 1.95),
                              min_size = 3,
                              distance_threshold = 1,
                              mu_sd = c(0.25, 0.5, 1, 2),
                              noise_sd = c(1),
                              bandwidth_sd = c(1, 3),
                              tyk_exp = c(1, 0.5))
configurations <- expand.grid(seed = 1:100,
                              dims = list(c(200, 1), c(14, 14)),
                              threshold = c(1.65, 1.95),
                              min_size = 3,
                              distance_threshold = 1,
                              mu_sd = c(0.25, 0.5, 1, 2),
                              noise_sd = c(1),
                              bandwidth_sd = c(1, 3),
                              tyk_exp = c(2, 1.5, 0.25))

registerDoParallel(4)
nothing <- foreach(i = 1:nrow(configurations)) %dopar% {
  result <- NULL
  try(result <- run_sim(configurations[i, ]))
  if(is.null(result)) {
    a <- c()
    filename <- sprintf("simulations/001/errors/iter_%s_%s_failed.rds", i, "B")
    saveRDS(a, file = filename)
  } else {
    filename <- sprintf("simulations/001/results/iter_%s_%s.rds", i, "B")
    saveRDS(result, file = filename)
  }
  return(result)
}





