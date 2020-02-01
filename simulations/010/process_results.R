library(data.table)
library(ggplot2)
library(magrittr)

sim_results_dir <- "simulations/010/results"
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

# Coverage rate -------------
cover_cols <- c("sim_id", "mu_sd", "bandwidth_sd", "threshold", "rho", "n_selected", "naive_cover", "profile_cover", "dims", "setting")
coverage_rate <- results[, .SD, .SDcols = cover_cols]
coverage_rate <- melt(coverage_rate, id.vars = setdiff(cover_cols, c("profile_cover", "naive_cover")),
                      variable.name = "ci_method", value.name = "coverage")
coverage_rate <- coverage_rate[, .(cover = mean(coverage),
                                   wcover = weighted.mean(coverage, n_selected)),
                         by = .(sim_id, mu_sd, bandwidth_sd, threshold, rho, ci_method, dims, setting)]
coverage_rate <- coverage_rate[, .(cover = mean(cover), coverSD = sd(cover) / sqrt(.N),
                                   wcover = mean(wcover), wcoverSD = sd(wcover) / sqrt(.N)),
                               by = .(mu_sd, bandwidth_sd, threshold, rho, ci_method, dims, setting)]
coverage_rate <- rbind(coverage_rate[, .(mu_sd, bandwidth_sd, setting, threshold, rho, ci_method, dims, cover = cover, coverSD = coverSD, avg = "simple")],
                       coverage_rate[, .(mu_sd, bandwidth_sd, setting, threshold, rho, ci_method, dims, cover = wcover, coverSD = wcoverSD, avg = "weighted")])
ci_quant <- qnorm(1 - 0.05 / nrow(coverage_rate))
ggplot(coverage_rate[ci_method == "profile_cover"],
       aes(x = mu_sd, y = cover, col = setting, linetype = avg)) +
  geom_line() +
  theme_bw() +
  facet_grid(dims + threshold ~ rho + bandwidth_sd, scales = "free", labeller = "label_both") +
  geom_hline(yintercept = 0.95) +
  geom_segment(aes(xend = mu_sd, y = cover - ci_quant * coverSD,
                   yend = cover + ci_quant * coverSD))
ggsave("simulations/010/figures/ci_coverage_rate.pdf")

# What are our misses? ------------------
misses <- results[, .(sim_id, mu_sd, bandwidth_sd, threshold, rho, dims, true_mean,
                      upper = abs(true_mean) > profile_uci * sign(true_mean) & !profile_cover,
                      lower = abs(true_mean) < profile_uci * sign(true_mean) & !profile_cover,
                      overall = 1 - profile_cover)]
misses <- misses[, .(upper = mean(upper),
                     lower = mean(lower),
                     overall = mean(overall)),
                 by = .(sim_id, mu_sd, bandwidth_sd, threshold, rho, dims)]
misses <- misses[, .(upper = mean(upper),
                     lower = mean(lower),
                     overall = mean(overall)),
                 by = .(mu_sd, bandwidth_sd, threshold, rho, dims)]
misses <- melt(misses, id.vars = setdiff(names(misses), c("lower", "upper", "overall")),
               variable.name = "direction",
               value.name = "percentage")
ggplot(misses, aes(x = mu_sd, y = percentage, col = direction, linetype = direction)) +
  facet_grid(dims + threshold ~ rho + bandwidth_sd, scales = "free", labeller = "label_both") +
  theme_bw() +
  geom_line() +
  geom_hline(yintercept = 0.05) +
  geom_hline(yintercept = 0.025, linetype = 2)

# Power ---------
power <- results[, .(sim_id, mu_sd, bandwidth_sd, threshold, rho, dims, setting,
                     naive_power = 0 < naive_lci | 0 > naive_uci,
                     profile_power = 0 < profile_lci | 0 > profile_uci)]
power <- power[, .(naive_power = mean(naive_power),
                   profile_power = mean(profile_power)),
               by = .(sim_id, mu_sd, bandwidth_sd, threshold, rho, dims, setting)]
power <- power[, .(naive_power = mean(naive_power),
                   profile_power = mean(profile_power)),
               by = .(mu_sd, bandwidth_sd, threshold, rho, dims, setting)]
power <- melt(power, id.vars = setdiff(names(power), c("naive_power", "profile_power")),
              variable.name = "method", value.name = "power")
ggplot(power, aes(x = mu_sd, y = power, linetype = method, col = setting)) +
  geom_line() +
  theme_bw() +
  facet_grid(dims + threshold ~ rho + bandwidth_sd, scales = "free", labeller = "label_both") +
  ylim(0, 1)









