library(dplyr)
library(reshape2)
library(ggplot2)
library(selectiveROI)
library(data.table)

writelog <- function(str, ...) {
  str = sprintf(str, ...)
  str = sprintf("%s: %s", Sys.time(), str)
  print(str)
}

generate_experiment_data <- function(p, rho = 0.9,
                                     mu.sd = 4,
                                     smooth.df = 0.2,
                                     threshold = 1.5,
                                     seed = NULL) {
  # Generating Data ---------------
  if (!is.null(seed)) {
    set.seed(seed)
  }
  sigma <- rho^as.matrix(dist(1:p))
  mu <- predict(smooth.spline(rnorm(p, sd = mu.sd), df = p * smooth.df))$y
  y <- as.numeric(mvtnorm::rmvnorm(1, mean = mu, sigma = sigma))
  threshold <- 1.5 / sqrt(diag(sigma))
  select <- abs(y) > threshold

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
  return(list(y = y,
              mu = mu,
              sigma = sigma,
              select = select,
              threshold = threshold))
}





# sampling -----
# probability under true param
# Computing MLE using old MLE function -----
raw.res = list()
summ.res = list()
experiments = expand.grid(p = 100,
                          rho = c(0.6, 0.9),
                          mu.sd = c(2, 4, 6),
                          smooth.df = c(0.2, 0.4),
                          threshold = c(1.64, 1.96),
                          seed = 1:100)
slacks = c(0, 0.1, 0.2, 0.5, 1, 2)
for (j in 1:nrow(experiments)) {
  params = lapply(experiments[j, ], function(x) x)
  dat.ok = FALSE
  set.seed(experiments$seed[j])
  params$seed = NULL
  while(!dat.ok) {
    dat = do.call(generate_experiment_data, params)
    list2env(dat, envir = parent.frame())
    dat.ok = sum(select) >= 3
    dat.ok = dat.ok & all(sign(y[select]) == sign(y[select][1]))
  }
  res = list()
  mmu = mean(mu[select])
  y =  y * sign(mmu)
  mu = mu * sign(mmu)
  res[[paste0("iter_", j, "_naive")]] = y
  res[[paste0("iter_", j, "_true")]] = mu
  res[[paste0("iter_", j, "_select")]] = select
  for (i in 1:length(slacks)) {
    writelog("iter %s, slack %s, n-selected %s out of %s in cluster",
             j, slacks[i], sum(select), length(select))
    mleest = NULL
    try(mleest <- roiMLE(y, sigma, threshold, compute = "mle",
                     coordinates = matrix(1:length(y), ncol = 1),
                     selected = select,
                     regularization_slack = slacks[i]))
    medest = NULL
    try(medest <- roiMLE(y, sigma, threshold, compute = "upper-CI",
                     ci_alpha = 0.49,
                     coordinates = matrix(1:length(y), ncol = 1),
                     selected = select,
                     regularization_slack = slacks[i]))
    res[[paste0("iter_", j, "_slack_", slacks[i], "_mle")]] = mleest$conditional
    res[[paste0("iter_", j, "_slack_", slacks[i], "_med")]] = medest$conditional
  }

  bias <- sapply(res, function(x, m, s) mean(x[s]) - mean(m[s]),
                 mu, select)
  avg_mse <- sapply(res, function(x, m, s) (mean(x[s]) - mean(m[s]))^2,
                    mu, select)
  vec_mse <- sapply(res, function(x, m, s) mean((x[s] - m[s])^2),
                    mu, select)
  iter_names <- names(res)
  tab <- data.table(name = iter_names,
                    avg_mse = avg_mse,
                    vec_mse = vec_mse,
                    bias = bias)
  summ.res[[j]] = tab
  raw.res[[j]] = res
  print(tab)
}
# saveRDS(list(summ = summ.res, raw = raw.res),
#         "simulations/011/list_of_all_B.rds")

tab = rbindlist(summ.res)
tab[, estimate := "naive"]
tab[grepl("mle", name), estimate := "mle"]
tab[grepl("med", name), estimate := "median"]
tab[grepl("true", name), estimate := "true"]
tab[grepl("select", name), estimate := "select"]
tab = tab[!(estimate %in% c("true", "select"))]
tab.slack = tab[estimate != "naive"]$name %>% strsplit("slack_") %>%
  sapply(function(x) x[[2]]) %>% strsplit("_") %>% sapply(function(x) x[1]) %>%
  as.numeric()
tab[estimate != "naive", slack := tab.slack]
tab[estimate == "naive", slack := 1]
tab = melt(tab, id.vars = c("name", "estimate", "slack"),
           variable.name = "stat")
tab.iters = tab$name %>% strsplit("iter_") %>% sapply(function(x) x[[2]]) %>%
  strsplit("_") %>% sapply(function(x) x[[1]]) %>% as.numeric()
tab[, iter := tab.iters]
setDT(experiments)
experiments[, iter := .I]
tab = merge(tab, experiments, by = c("iter"))
tab[, iter := NULL]
bycols = setdiff(names(tab), c("value", "name", "seed"))
by.expr = parse(text = sprintf("list(%s)", paste0(bycols, collapse = ", ")))
summ = tab[, list(value = mean(unlist(value))),
           by = eval(by.expr)]
summ = summ[is.finite(value)]
summ[, smooth.df := smooth.df]

ggplot(summ[estimate %in% c("mle")],
       aes(x = slack, y = value, col = factor(mu.sd), linetype = estimate)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  facet_wrap(~ stat, ncol = 2) +
  facet_grid(stat ~ smooth.df + rho + threshold, scales = "free", labeller = "label_both") +
  geom_hline(yintercept = 0, linetype = 2)


# Raw ----
raw = lapply(1:length(raw.res), function(i) {
  x = raw.res[[i]]
  if(is.null(x)) return(NULL)
  names = names(x)
  x = as.data.table(do.call("rbind", x))
  x[, name := names]
  x = melt(x, id.vars = "name")
})
raw = rbindlist(raw)
iters = raw$name %>% strsplit("_") %>% sapply(function(x) x[[2]]) %>% as.numeric()
raw[, iter := iters]
ests = strsplit(raw$name, "_") %>% sapply(tail, 1)
raw[, estimate := ests]
not.est = raw[estimate %in% c("naive", "true", "select")]
not.est = dcast(not.est, variable + iter ~ estimate)
raw = raw[estimate %in% c("median", "mle")]
raw = merge(raw, not.est, by = c("variable", "iter"))
raw = raw[select == 1]
summ = raw[, .(est = mean(value), naive = mean(naive), true = mean(true)),
           by = .(name)]
summ[sample(nrow(summ))][1:10]
