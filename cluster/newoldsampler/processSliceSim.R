filenames <- as.list(dir(path = 'cluster/newoldsampler/results/', pattern="sliceComparison_C_*")) # W Dependence
filenames <- lapply(filenames, function(x) paste0('cluster/newoldsampler/results/', x))
filenames <- as.character(filenames)
results <- vector(length(filenames), mode = "list")
for(i in 1:length(filenames)) {
  results[[i]] <- readRDS(filenames[i])
}

results <- do.call("c", results)
results <- do.call("c", results)
results <- results[sapply(results, function(x) !is.null(x))]

# runtime ------
computeRuntime <- function(x) {
  dat <- c(rho = x$params[[2]], size = length(x$mu),
           gibbs = x$times[[1]],
           slice = x$times[[2]])
}
runtime <- sapply(results, computeRuntime) %>% t()
runtime <- data.frame(runtime)

# runtime <- group_by(runtime, size) %>%
#   summarize(gibbs = median(gibbs), slice = median(slice))
runtime <- melt(runtime, id = c("size", "rho"))
names(runtime)[3:4] <- c("method", "time")
ggplot(runtime) +
  geom_smooth(aes(x = size, y = time, col = method), method = "loess") +
  theme_bw() +
  ylab("time in secs") +
  xlab("cluster size")

# estimation error -----
computeError <- function(x) {
  select <- x$selected
  signs <- sign(x$observed[x$selected][1])
  dat <- c(musd = x$params[[1]],
           rho = x$params[[2]],
           threshold = x$params[[3]],
           size = sum(x$selected),
           observed = signs * (mean(x$observed[select]) - mean(x$mu[select])),
           gibbs = signs * (mean(x$oldmle[select]) - mean(x$mu[select])),
           slice = signs * (mean(x$newmle[select]) - mean(x$mu[select])))
}
error <- sapply(results, computeError) %>% t()
error <- data.frame(error)
error <- melt(error, c("musd", "rho", "threshold", "size"))
bins <- c(1, 5, 10, 15, 20, 30, Inf)
binnames <- as.character(bins)
error$bsize <- length(bins) + 1  -sapply(error$size, function(x) sum(x <= bins))
error$sizeFactor <- factor(binnames[error$bsize], levels = binnames)
names(error)[5:6] <- c("method", "error")
ggplot(subset(error, rho == 0.9), aes(x = sizeFactor, y = abs(error), col = method)) +
  geom_boxplot() + #geom_jitter(size = 0.1) +
  theme_bw() +
  facet_grid(musd ~ threshold, labeller = "label_both") +
  geom_hline(yintercept = 0) +
  xlab("< size <=")

ggplot(subset(error, rho == 0.9), aes(x = sizeFactor, y = error, col = method)) +
  geom_boxplot() + #geom_jitter(size = 0.1) +
  theme_bw() +
  facet_grid(musd ~ threshold, labeller = "label_both") +
  geom_hline(yintercept = 0) +
  xlab("< size <=")

ggplot(error, aes(x = factor(musd), y = error, col = method)) +
  geom_boxplot() + #geom_jitter(size = 0.1) +
  theme_bw() +
  facet_grid(rho ~ threshold, labeller = "label_both") +
  geom_hline(yintercept = 0)

ggplot(error, aes(x = factor(musd), y = abs(error), col = method)) +
  geom_boxplot() + #geom_jitter(size = 0.1) +
  theme_bw() +
  facet_grid(rho ~ threshold, labeller = "label_both") +
  geom_hline(yintercept = 0)

