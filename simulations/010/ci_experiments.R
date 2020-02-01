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

      # Computing MLE --------------------
      # mle <- roiMLE(y = y[this_cluster],
      #               cov = this_covariance, threshold,
      #               coordinates = coordinates[this_cluster, , drop = FALSE],
      #               selected = selected,
      #               regularization_slack = 1,
      #               progress = TRUE,
      #               mle_control = control)

      # Computing ci --------------------
      control <- roi_mle_control(grad_iterations = grad_iterations,
                                 step_size_coef = grad_step_size,
                                 step_rate = grad_step_rate,
                                 samp_per_iter = config[["samp_per_iter"]],
                                 grad_delay = delay,
                                 assume_convergence = assume_convergence,
                                 impute_boundary = config[["imputation_method"]],
                                 RB_mult = config[["RB_mult"]])
      if(mean(y[this_cluster][selected]) > 0) {
        control$grad_iterations <- 4500
      } else {
        control$grad_iterations <- 1000
      }
      lower_ci_fit <- roiMLE(y = y[this_cluster],
                             cov = this_covariance, threshold,
                             compute = "lower-CI",
                             ci_alpha = 0.025,
                             coordinates = coordinates[this_cluster, , drop = FALSE],
                             selected = selected,
                             regularization_slack = 1,
                             progress = TRUE,
                             mle_control = control)
      n_steps <- length(lower_ci_fit$ci_lim_path)
      ci_comp_index <- floor(n_steps * 2/3):n_steps
      lower_ci <- mean(lower_ci_fit$ci_lim_path[ci_comp_index])

      if(mean(y[this_cluster][selected]) < 0) {
        control$grad_iterations <- 4500
      } else {
        control$grad_iterations <- 1500
      }
      upper_ci_fit <- roiMLE(y = y[this_cluster],
                             cov = this_covariance, threshold,
                             compute = "upper-CI",
                             ci_alpha = 0.025,
                             coordinates = coordinates[this_cluster, , drop = FALSE],
                             selected = selected,
                             regularization_slack = 1,
                             progress = TRUE,
                             mle_control = control)
      n_steps <- length(upper_ci_fit$ci_lim_path)
      ci_comp_index <- floor(n_steps * 2/3):n_steps
      upper_ci <- mean(upper_ci_fit$ci_lim_path[ci_comp_index])


      # Reporting
      naive_cov <- covariance[clusters[[i]], clusters[[i]]]
      obs_mean <- mean(y[this_cluster][selected])
      avg_vec <- rep(1, sum(selected)) / sum(selected)
      naive_sd <- sqrt(as.numeric(t(avg_vec) %*% naive_cov %*% avg_vec))
      naive_ci <- c(obs_mean - 1.96 * naive_sd, obs_mean + 1.96 * naive_sd)
      true_mean <- true_mean <- mean(mu[this_cluster][selected])
      iter_results[[i]] <- c(n_selected = sum(selected),
                             size = length(this_cluster),
                             true_mean = true_mean,
                             obs_mean = obs_mean,
                             naive_lci = naive_ci[1],
                             naive_uci = naive_ci[2],
                             naive_cover = naive_ci[1] <= true_mean & naive_ci[2] >= true_mean,
                             profile_lci = lower_ci,
                             profile_uci = upper_ci,
                             profile_cover = lower_ci <= true_mean & upper_ci >= true_mean)
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
                                grad_step_size = c(0.5),
                                grad_iterations = c(2100),
                                samp_size = c(10^4),
                                tyk_exp = c(1),
                                min_size = 3,
                                distance_threshold = 1,
                                n_chains = c(5),
                                burnin = c(1000),
                                mu_sd = c(0.25, 0.5),
                                bandwidth_sd = c(1, 3),
                                seed = 1:100,
                                threshold = c(1.65, 1.95),
                                rho = c(0.8, 0.6),
                                grad_step_rate = c(0.55),
                                noise_sd = c(1),
                                dims = c("1D", "2D"),
                                samp_per_iter = c(100),
                                RB_mult = 1,
                                setting = "baseline")
configurations_B <- expand.grid(imputation_method = c("smooth"),
                                grad_step_size = c(0.5),
                                grad_iterations = c(2100),
                                samp_size = c(10^4),
                                tyk_exp = c(1),
                                min_size = 3,
                                distance_threshold = 1,
                                n_chains = c(5),
                                burnin = c(1000),
                                mu_sd = c(0.25, 0.5),
                                bandwidth_sd = c(1, 3),
                                seed = 1:100,
                                threshold = c(1.65, 1.95),
                                rho = c(0.8, 0.6),
                                grad_step_rate = c(0.55),
                                noise_sd = c(1),
                                RB_mult = c(0.5),
                                dims = c("1D", "2D"),
                                samp_per_iter = c(100),
                                setting = "smaller_RB")
configurations_C <- expand.grid(imputation_method = c("smooth"),
                                grad_step_size = c(0.5),
                                grad_iterations = c(2100),
                                samp_size = c(10^4),
                                tyk_exp = c(1),
                                min_size = 3,
                                distance_threshold = 1,
                                n_chains = c(5),
                                burnin = c(1000),
                                mu_sd = c(0.25, 0.5),
                                bandwidth_sd = c(1, 3),
                                seed = 1:100,
                                threshold = c(1.65, 1.95),
                                rho = c(0.8, 0.6),
                                grad_step_rate = c(0.55),
                                noise_sd = c(1),
                                RB_mult = c(2),
                                dims = c("1D", "2D"),
                                samp_per_iter = c(100),
                                setting = "largerRB")
configurations_D <- expand.grid(imputation_method = c("smooth"),
                                grad_step_size = c(0.5),
                                grad_iterations = c(2100),
                                samp_size = c(10^4),
                                tyk_exp = c(1),
                                min_size = 3,
                                distance_threshold = 1,
                                n_chains = c(5),
                                burnin = c(1000),
                                mu_sd = c(0.25, 0.5),
                                bandwidth_sd = c(1, 3),
                                seed = 1:100,
                                threshold = c(1.65, 1.95),
                                rho = c(0.8, 0.6),
                                grad_step_rate = c(0.55),
                                noise_sd = c(1),
                                RB_mult = c(1),
                                dims = c("1D", "2D"),
                                samp_per_iter = c(200),
                                setting = "more_samp")
configurations_E <- expand.grid(imputation_method = c("smooth"),
                                grad_step_size = c(0.5),
                                grad_iterations = c(2100),
                                samp_size = c(10^4),
                                tyk_exp = c(1),
                                min_size = 3,
                                distance_threshold = 1,
                                n_chains = c(5),
                                burnin = c(1000),
                                mu_sd = c(0.25, 0.5),
                                bandwidth_sd = c(1, 3),
                                seed = 1:100,
                                threshold = c(1.65, 1.95),
                                rho = c(0.8, 0.6),
                                grad_step_rate = c(0.85),
                                noise_sd = c(1),
                                RB_mult = c(1),
                                dims = c("1D", "2D"),
                                samp_per_iter = c(100),
                                setting = "smaller_step_rate")
configurations_F <- expand.grid(imputation_method = c("smooth"),
                                grad_step_size = c(0.5),
                                grad_iterations = c(4200),
                                samp_size = c(10^4),
                                tyk_exp = c(1),
                                min_size = 3,
                                distance_threshold = 1,
                                n_chains = c(5),
                                burnin = c(1000),
                                mu_sd = c(0.25, 0.5),
                                bandwidth_sd = c(1, 3),
                                seed = 1:100,
                                threshold = c(1.65, 1.95),
                                rho = c(0.8, 0.6),
                                grad_step_rate = c(0.55),
                                noise_sd = c(1),
                                RB_mult = c(1),
                                dims = c("1D", "2D"),
                                samp_per_iter = c(100),
                                setting = "more_iter")
configurations_G <- expand.grid(imputation_method = c("smooth"),
                                grad_step_size = c(0.5),
                                grad_iterations = c(3150),
                                samp_size = c(10^4),
                                tyk_exp = c(1),
                                min_size = 3,
                                distance_threshold = 1,
                                n_chains = c(5),
                                burnin = c(1000),
                                mu_sd = c(0.25, 0.5),
                                bandwidth_sd = c(1, 3),
                                seed = 1:100,
                                threshold = c(1.65, 1.95),
                                rho = c(0.8, 0.6),
                                grad_step_rate = c(0.85),
                                noise_sd = c(1),
                                RB_mult = c(4),
                                dims = c("1D", "2D"),
                                samp_per_iter = c(400),
                                setting = "+samp++RB-step+iter")

configurations <- rbind(configurations_A,
                        configurations_B,
                        configurations_C,
                        configurations_D,
                        configurations_E,
                        configurations_F,
                        configurations_G)
# set.seed(1)
# configurations = configurations[sample(1:nrow(configurations)), ]
configurations$dims <- factor(configurations$dims, levels = c("1D", "2D"))
registerDoParallel(4)
sim_path <- "simulations/010"
have_results <- dir(sprintf("%s/results", sim_path)) %>% strsplit("_") %>%
  sapply(function(x) x[[2]]) %>% as.numeric()
to_run <- setdiff(1:nrow(configurations), have_results)
foreach(i = to_run) %dopar% {
  log_file <- sprintf("%s/errors/iter_%s_%s_log.txt", sim_path, i, "C")
  error_file <- sprintf("%s/errors/iter_%s_%s_error.txt", sim_path, i, "C")
  sink(log_file)
  result <- NULL
  try(result <- run_sim(configurations[i, ]))
  if(!is.null(result)) {
    file.remove(log_file)
    filename <- sprintf("%s/results/iter_%s_%s.rds", sim_path, i, "C")
    saveRDS(result, file = filename)
    sink()
  } else {
    file.rename(log_file, error_file)
    sink()
  }
  return(NULL)
}





