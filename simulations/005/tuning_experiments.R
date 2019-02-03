library(foreach)
library(magrittr)
library(doParallel)
library(selectiveROI)
library(data.table)

# simulation parameters ---------
run_sim <- function(config) {
  set.seed(config[["seed"]])
  if(config[["dims"]] == "1D") {
    dims <- c(200, 1)
  } else if(config[["dims"]] == "2D") {
    dims <- c(14, 14)
  }

  data <- generate_smooth_coordinate_data(n = 1,
                                          dims = dims,
                                          rho = config[["rho"]],
                                          noise_sd = config[["noise_sd"]],
                                          mu_sd = config[["mu_sd"]],
                                          bandwidth_sd = config[["bandwidth_sd"]],
                                          dist_method = "manhattan")
  grad_step_size <- config[["grad_step_size"]]
  grad_step_rate <- config[["grad_step_rate"]]
  grad_iterations <- config[["grad_iterations"]]
  assume_convergence <- floor(grad_iterations / 3)
  delay <- floor(assume_convergence / 3)
  samp_size <- config[["samp_size"]]
  threshold <- config[["threshold"]]
  lims <- max(abs(c(data$y[1, ], data$mu[1, ])))
  lims <- c(-lims, lims)

  if(FALSE) {
    plot(data$y[1, ], type = "l", ylim = lims)
    lines(data$mu[1, ], type = "l", col = "red")
    abline(h = c(threshold, - threshold), col = "grey", lty = 2)
  }

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
      if(config[["imputation_method"]] == "none") {
        boundary <- c()
      } else {
        boundary <- setdiff(1:n_coordinates, this_cluster)[boundary]
      }
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
                    stepSizeCoef = grad_step_size,
                    stepRate = grad_step_rate,
                    sampPerIter = config[["samp_per_iter"]],
                    delay = delay,
                    maxiter = grad_iterations,
                    assumeConvergence = assume_convergence,
                    nsamp = 0,
                    imputeBoundary = config[["imputation_method"]],
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
                            stepSizeCoef = grad_step_size,
                            stepRate = grad_step_rate,
                            sampPerIter = config[["samp_per_iter"]],
                            delay = delay,
                            maxiter = grad_iterations,
                            assumeConvergence = assume_convergence,
                            nsamp = samp_size,
                            imputeBoundary = config[["imputation_method"]],
                            progress = TRUE)
      samp <- rowMeans(profile_mle$sample[, selected])
      quants <- quantile(samp, probs = c(0.025, 0.975))
      covered <- quants[1] < obs_mean & quants[2] > obs_mean
      naive_cov <- covariance[clusters[[i]], clusters[[i]]]
      avg_vec <- rep(1, sum(selected)) / sum(selected)
      naive_sd <- sqrt(as.numeric(t(avg_vec) %*% naive_cov %*% avg_vec))
      naive_ci <- c(obs_mean - 1.96 * naive_sd, obs_mean + 1.96 * naive_sd)
      naive_cover <- true_mean >= naive_ci[1] & true_mean <= naive_ci[2]
      iter_results[[i]] <- c(n_selected = sum(selected),
                             size = length(this_cluster),
                             true_mean = true_mean,
                             obs_mean = obs_mean,
                             mle = mean(mle$conditional[selected]),
                             q025 = quants[1],
                             q975 = quants[2],
                             covered = covered,
                             naive_cover = naive_cover,
                             imp_cover = naive_cover | covered)
      print(iter_results[[i]])
      a <- 0
    }
    result_list[[m]] <- do.call("rbind", iter_results)
    lapply(result_list[1:m], colMeans) %>% do.call("rbind", .) %>% colMeans() %>% print()
  }
  return(list(parameters = config,
              results = result_list))
}

configurations_A <- expand.grid(imputation_method = c("smooth"),
                                samp_per_iter = c(20, 50, 100, 200),
                                grad_step_size = c(0.5),
                                grad_iterations = c(2100),
                                samp_size = c(10^4 * 2),
                                tyk_exp = c(1),
                                min_size = 3,
                                distance_threshold = 1,
                                mu_sd = c(0.25, 0.5, 1, 2),
                                bandwidth_sd = c(1, 3),
                                seed = 1:100,
                                threshold = c(1.65, 1.95),
                                rho = c(0.9, 0.6),
                                grad_step_rate = c(0.55, 0.75),
                                noise_sd = c(1),
                                dims = c("1D", "2D"))
configurations <- configurations_A
configurations$dims <- factor(configurations$dims, levels = c("1D", "2D"))
registerDoParallel(4)
sim_path <- "simulations/005"
have_results <- dir(sprintf("%s/results", sim_path)) %>% strsplit("_") %>%
  sapply(function(x) x[[2]]) %>% as.numeric()
to_run <- setdiff(1:nrow(configurations), have_results)
foreach(i = to_run) %dopar% {
  log_file <- sprintf("%s/errors/iter_%s_%s_log.txt", sim_path, i, "A")
  error_file <- sprintf("%s/errors/iter_%s_%s_error.txt", sim_path, i, "A")
  sink(log_file)
  result <- NULL
  try(result <- run_sim(configurations[i, ]))
  if(!is.null(result)) {
    file.remove(log_file)
    filename <- sprintf("%s/results/iter_%s_%s.rds", sim_path, i, "A")
    saveRDS(result, file = filename)
    sink()
  } else {
    file.rename(log_file, error_file)
    sink()
  }
  return(NULL)
}





