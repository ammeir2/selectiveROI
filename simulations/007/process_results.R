library(data.table)
library(ggplot2)
library(magrittr)

sim_results_dir <- "simulations/007/results"
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
                                      "imp_cover.2.5%", "size", "n_chains"))
estimation_error <- results[, .SD, .SDcols = est_cols]
estimation_error <- melt(estimation_error,
                         id.vars = setdiff(est_cols, c("obs_mean", "mle")),
                         value.name = "estimate", variable.name = "method")
estimation_error <- estimation_error[, .(rmse = sqrt(mean((estimate - true_mean)^2)),
                                         wrmse = sqrt(weighted.mean((estimate - true_mean)^2,
                                                                    w = n_selected))),
                                     by = .(sim_id, imputation_method, grad_step_rate, grad_step_size,
                                            grad_iterations, threshold, mu_sd, noise_sd,
                                            bandwidth_sd, rho, dims, method, samp_per_iter)]
estimation_error <- melt(estimation_error,
                         id.vars = setdiff(names(estimation_error), c("rmse", "wrmse")),
                         variable.name = "error_type", value.name = "error")
estimation_error <- estimation_error[, .(mean_error = mean(error), sd_error = sd(error) / sqrt(.N)),
                                     by = .(imputation_method, grad_step_rate, grad_step_size,
                                            grad_iterations, threshold, mu_sd, noise_sd,
                                            bandwidth_sd, rho, dims, method, error_type, samp_per_iter)]
ggplot(estimation_error, aes(x = mu_sd, y = mean_error, col = method, linetype = error_type)) +
  facet_grid(dims + threshold ~ rho + bandwidth_sd, scales = "free", labeller = "label_both") +
  theme_bw() +  geom_line()
ggsave("simulations/007/figures/estimation_error.pdf")

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
                                threshold, mu_sd, bandwidth_sd, rho, dims, ci_method, samp_per_iter,
                                n_chains, burnin)]
coverage_rate <- coverage_rate[, .(cover = mean(cover), coverSD = sd(cover) / sqrt(.N),
                                   wcover = mean(wcover), wcoverSD = sd(wcover) / sqrt(.N)),
                               by = .(grad_step_size, grad_step_rate, grad_iterations, samp_size, tyk_exp,
                                      threshold, mu_sd, bandwidth_sd, rho, dims, ci_method, samp_per_iter,
                                      n_chains, burnin)]
# coverage_rate <- melt(coverage_rate,
#                       id.vars = setdiff(names(coverage_rate), c("cover", "wcover")),
#                       variable.name = "measure", value.name = "coverage")
coverage_rate[, tune := sprintf("rate%s_gardsamps%s_burnin%s_nchains%s",
                                grad_step_rate, samp_per_iter, burnin, n_chains)]
coverage_rate[, coverage := cover]
coverage_rate[, coverageSD := coverSD]
ggplot(coverage_rate[ci_method == "profile" & burnin == 1000],
       aes(x = mu_sd, y = coverage, col = tune, linetype = tune)) +
  geom_line() +
  theme_bw() +
  facet_grid(dims + threshold ~ rho + bandwidth_sd, scales = "free", labeller = "label_both") +
  geom_hline(yintercept = 0.95) +
  geom_segment(aes(xend = mu_sd, y = coverage - 2 * coverageSD,
                   yend = coverage + 2 * coverageSD))
ggsave("simulations/006/figures/ci_coverage_rate.pdf")


















