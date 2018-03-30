filenames <- as.list(dir(path = 'cluster/cubes/results/', pattern="cubeROI_A_*"))
filenames <- lapply(filenames, function(x) paste0('cluster/cubes/results/', x))
filenames <- as.character(filenames)
results <- vector(length(filenames), mode = "list")
for(i in 1:length(filenames)) {
  results[[i]] <- readRDS(filenames[i])
}
results <- do.call("c", results)
results <- do.call("c", results)

# Type I Error -------
compTypeI <- function(res) {
  nclusts <- length(res)
  # print(nclusts)
  results <- NULL
  try(results <- data.frame(id = runif(1),
                            clust = 1:length(res),
                            size = sapply(res, function(x) x[[1]]["size"]),
                            rho = sapply(res, function(x) x[[1]]["rho"]),
                            threshold = sapply(res, function(x) x[[1]]["BHlevel"]),
                            reject = sapply(res, function(x) x[[2]][3] <= 0.025)))
  if(is.null(results)) {
    return(NULL)
  }

  return(results)
}
zerosnr <- results[sapply(results, function(x) x[[1]][[1]]["snr"] == 0)]
typeI <- lapply(zerosnr, compTypeI)
typeI <- data.table::rbindlist(typeI)
typeI <- data.frame(typeI)

sumTypeI <- group_by(typeI, id, rho, threshold) %>%
  summarize(typeI = weighted.mean(reject, size))
sumTypeI <- group_by(sumTypeI, rho, threshold) %>%
  summarize(typeIsd = sd(typeI) / sqrt(length(typeI)),
            typeI = mean(typeI)) %>% as.data.frame()
sumTypeI

# Coverage rate ---------
library(magrittr)
library(dplyr)
library(ggplot2)
compCover <- function(res) {
  nclusts <- length(res)
  # print(nclusts)
  results <- NULL
  try(results <- data.frame(id = runif(1),
                        clust = 1:length(res),
                        size = sapply(res, function(x) x[[1]]["size"]),
                        snr = sapply(res, function(x) x[[1]]["snr"]),
                        rho = sapply(res, function(x) x[[1]]["rho"]),
                        threshold = sapply(res, function(x) x[[1]]["BHlevel"]),
                        cover = sapply(res, function(x) x[[2]][2] >= 0.05)))
  if(is.null(results)) {
    return(NULL)
  }

  return(results)
}
cover <- sapply(results, compCover)
cover <- data.table::rbindlist(cover)
cover <- data.frame(cover)

sumcover <- group_by(cover, id, rho, snr, threshold) %>%
  summarize(cover = weighted.mean(cover, size))
sumcover <- group_by(sumcover, rho, threshold, snr) %>%
  summarize(coversd = sd(cover) / sqrt(length(cover)),
            cover = mean(cover)) %>% as.data.frame()

ggplot(sumcover, aes(x = snr, y = cover, col = factor(rho))) +
  geom_line() +
  geom_segment(aes(xend = snr, y = cover - 2*coversd, yend = pmin(cover + 2*coversd, 1))) +
  facet_wrap(~ threshold) +
  geom_hline(yintercept = 0.95) +
  theme_bw()

