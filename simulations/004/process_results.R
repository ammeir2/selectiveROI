library(data.table)
library(ggplot2)
library(magrittr)

sim_results_dir <- "simulations/004/results"
files <- dir(sim_results_dir)
files <- paste(sim_results_dir, "/", files, sep = "")
param_list <- vector(length(files), mode = "list")
result_list <- vector(length(files), mode = "list")
print(sprintf("processing %s files", length(files)))
for(i in 1:length(files)) {
  if(i %% 100 == 0) cat(i, " ")
  if(i %% 10^3 == 0) cat("\n")
  iter_result <- readRDS(files[i])
  result_tab <- as.data.table(iter_result$results[[1]])
  if(nrow(result_tab) == 0) next
  result_tab$sim_id <- i
  param_tab <- (iter_result$parameters)
  param_tab$sim_id <- i
  param_list[[i]] <- param_tab
  result_list[[i]] <- result_tab
}

results <- rbindlist(result_list)
parameters <- rbindlist(param_list)
results <- merge(parameters, results, by = "sim_id")
results[, seed := NULL]

# Estimation error ----------------
est_cols <- setdiff(names(results), c("samp_size", "tyk_exp", "min_size", "distance_threshold",
                                      "q025.2.5%", "q975.97.5%", "covered.2.5%", "naive_cover",
                                      "imp_cover.2.5%", "size"))
estimation_error <- results[, .SD, .SDcols = est_cols]
estimation_error <- melt(estimation_error,
                         id.vars = setdiff(est_cols, c("obs_mean", "mle")),
                         value.name = "estimate", variable.name = "method")
estimation_error <- estimation_error[, .(rmse = sqrt(mean((estimate - true_mean)^2)),
                                         wrmse = sqrt(weighted.mean((estimate - true_mean)^2,
                                                                    w = n_selected))),
                                     by = .(sim_id, imputation_method, grad_step_rate, grad_step_size,
                                            grad_iterations, threshold, mu_sd, noise_sd,
                                            bandwidth_sd, rho, dims, method)]
estimation_error <- melt(estimation_error,
                         id.vars = setdiff(names(estimation_error), c("rmse", "wrmse")),
                         variable.name = "error_type", value.name = "error")
estimation_error <- estimation_error[, .(mean_error = mean(error), sd_error = sd(error) / sqrt(.N)),
                                     by = .(imputation_method, grad_step_rate, grad_step_size,
                                            grad_iterations, threshold, mu_sd, noise_sd,
                                            bandwidth_sd, rho, dims, method, error_type)]
estimation_error[, tune := sprintf("rate%s_size%s_iter%s",
                                   grad_step_rate, grad_step_size, grad_iterations)]
ggplot(estimation_error, aes(x = mu_sd, y = mean_error, col = tune, linetype = method)) +
  facet_grid(dims + error_type ~ rho + bandwidth_sd + threshold, scales = "free", labeller = "label_both") +
  theme_bw() +  geom_line()

rankings <- estimation_error[order(dims, mu_sd, threshold, bandwidth_sd, rho, method, mean_error)][error_type == "rmse" & method == "mle"]
rankings <- rankings[, .(best = tune[1], worst = tune[.N]),
                     by = .(threshold, mu_sd, noise_sd, bandwidth_sd, rho, dims)]
rankings[order(best)]

# Coverage rate -------------
setnames(results, c("covered.2.5%", "naive_cover", "imp_cover.2.5%"), c("profile", "naive", "union"))
cover_cols <- setdiff(names(results), c("imputation_method", "min_size", "size", "true_mean",
                                        "obs_mean", "mle", "q025.2.5%", "q975.97.5%"))
coverage_rate <- results[, .SD, .SDcols = cover_cols]
coverage_rate <- melt(coverage_rate, id.vars = setdiff(cover_cols, c("profile", "naive", "union")),
                      variable.name = "ci_method", value.name = "coverage")
coverage_rate <- coverage_rate[, .(cover = mean(coverage),
                                   wcover = weighted.mean(coverage, n_selected)),
                         by = .(sim_id, grad_step_size, grad_step_rate, grad_iterations, samp_size, tyk_exp,
                                threshold, mu_sd, bandwidth_sd, rho, dims, ci_method)]
coverage_rate <- coverage_rate[, .(cover = mean(cover),
                                   wcover = mean(wcover)),
                               by = .(grad_step_size, grad_step_rate, grad_iterations, samp_size, tyk_exp,
                                      threshold, mu_sd, bandwidth_sd, rho, dims, ci_method)]
coverage_rate <- melt(coverage_rate,
                      id.vars = setdiff(names(coverage_rate), c("cover", "wcover")),
                      variable.name = "measure", value.name = "coverage")
coverage_rate[, tune := sprintf("size%s_rate%s_iter%s",
                                grad_step_size, grad_step_rate, grad_iterations)]
ggplot(coverage_rate[ci_method == "profile" & measure == "wcover"],
       aes(x = mu_sd, y = coverage, col = tune, linetype = factor(samp_size))) +
  geom_line() +
  theme_bw() +
  facet_grid(dims + threshold ~ rho + bandwidth_sd, scales = "free", labeller = "label_both") +
  geom_hline(yintercept = 0.95)


















