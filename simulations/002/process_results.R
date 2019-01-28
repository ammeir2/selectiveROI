library(data.table)
library(ggplot2)

sim_results_dir <- "simulations/002/results"
files <- dir(sim_results_dir)
files <- paste(sim_results_dir, "/", files, sep = "")
param_list <- vector(length(files), mode = "list")
result_list <- vector(length(files), mode = "list")
print(sprintf("processing %s files", length(files)))
for(i in 1:length(files)) {
  if(i %% 100) cat(i, " ")
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
results[, dim_type := "2D"]
results$dim_type[sapply(results$dims, function(x) x[1] == 200)] <- "1D"
results[, seed := NULL]
results[, dims := NULL]

# Estimation error ----------------
estimation_error <- results[, .(sim_id, dim_type, threshold, mu_sd, noise_sd,
                                imputation_method,
                                bandwidth_sd, n_selected,
                                obs_mean, mle, true_mean)]
estimation_error <- melt(estimation_error,
                         id.vars = c("sim_id", "threshold", "mu_sd", "noise_sd",
                                     "bandwidth_sd", "n_selected", "dim_type", "true_mean",
                                     "imputation_method"),
                         value.name = "estimate", variable.name = "method")
estimation_error <- estimation_error[, .(rmse = sqrt(mean((estimate - true_mean)^2)),
                                         wrmse = sqrt(weighted.mean((estimate - true_mean)^2,
                                                                    w = n_selected))),
                                     by = .(sim_id, threshold, mu_sd, noise_sd, bandwidth_sd,
                                            dim_type, method, imputation_method)]
estimation_error <- melt(estimation_error,
                         id.vars = c("sim_id", "threshold", "mu_sd", "noise_sd",
                                     "bandwidth_sd", "dim_type", "method",
                                     "imputation_method"),
                         variable.name = "error_type", value.name = "error")
estimation_error <- estimation_error[, .(mean_error = mean(error), sd_error = sd(error) / sqrt(.N)),
                                     by = .(threshold, mu_sd, noise_sd, bandwidth_sd, dim_type, method,
                                            error_type, imputation_method)]
ggplot(estimation_error, aes(x = mu_sd, y = mean_error, col = imputation_method, linetype = method)) +
  facet_grid(dim_type + error_type ~ bandwidth_sd + threshold, scales = "free", labeller = "label_both") +
  theme_bw() +
  geom_line()

# Coverage rate -------------
setnames(results, c("covered.2.5%", "naive_cover", "imp_cover.2.5%"), c("profile", "naive", "union"))
coverage_rate <- results[, .(size, sim_id, imputation_method, threshold, mu_sd,
                                   noise_sd, bandwidth_sd, profile, naive, union, dim_type)]
coverage_rate <- melt(coverage_rate, id.vars = c("size", "imputation_method", "threshold",
                                           "mu_sd", "noise_sd", "bandwidth_sd", "sim_id", "dim_type"),
                      variable.name = "ci_method", value.name = "coverage")
coverage_rate <- coverage_rate[, .(cover = mean(coverage),
                             wcover = weighted.mean(coverage, size)),
                         by = .(sim_id, threshold, mu_sd, noise_sd, bandwidth_sd, dim_type, imputation_method,
                                ci_method)]
coverage_rate <- coverage_rate[, .(cover = mean(cover),
                                   wcover = mean(wcover)),
                               by = .(threshold, mu_sd, noise_sd, bandwidth_sd, dim_type, imputation_method,
                                      ci_method)]
coverage_rate <- melt(coverage_rate,
                      id.vars = c("threshold", "mu_sd", "noise_sd", "bandwidth_sd", "dim_type", "imputation_method",
                                  "ci_method"),
                      variable.name = "measure", value.name = "coverage")
ggplot(subset(coverage_rate, ci_method == "profile"),
       aes(x = mu_sd, y = coverage, col = factor(imputation_method), linetype = ci_method)) +
  geom_line() +
  theme_bw() +
  facet_grid(dim_type + measure~ bandwidth_sd + threshold, scales = "free", labeller = "label_both") +
  geom_hline(yintercept = 0.95)


















