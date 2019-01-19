library(magrittr)
library(dplyr)
library(ggplot2)
library(reshape2)
# filenames <- as.list(dir(path = 'cluster/cubes/results/', pattern="cubeROI_E_*"))
# filenames <- as.list(dir(path = 'cluster/cubes/results/', pattern="cubeROI_F_*"))
# filenames <- as.list(dir(path = 'cluster/cubes/results/', pattern="cubeROI_G_*"))
# filenames <- as.list(dir(path = 'cluster/cubes/results/', pattern="cubeROI_H_*"))
# filenames <- as.list(dir(path = 'cluster/cubes/results/', pattern="cubeROI_I_*"))
# filenames <- as.list(dir(path = 'cluster/cubes/results/', pattern="cubeROI_J_*"))
# filenames <- as.list(dir(path = 'cluster/cubes/results/', pattern="flatROI_A_*"))
# filenames <- as.list(dir(path = 'cluster/cubes/results/', pattern="cubeROI_L_*"))
# filenames <- as.list(dir(path = 'cluster/cubes/results/', pattern="cubeROI_M_*"))
# filenames <- as.list(dir(path = 'cluster/cubes/results/', pattern="cubeROI_N_*"))
filenames <- as.list(dir(path = 'cluster/cubes/results/', pattern="cubeROI_O_*"))
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
  empty <- sapply(res, function(x) length(x) == 0)
  if(any(empty)) {
    res <- res[-which(empty)]
  }
  nclusts <- length(res)
  # print(nclusts)
  results <- NULL
  try(results <- data.frame(id = runif(1),
                            clust = 1:length(res),
                            size = sapply(res, function(x) x[[1]]["size"]),
                            rho = sapply(res, function(x) x[[1]]["rho"]),
                            threshold = sapply(res, function(x) x[[1]]["BHlevel"]),
                            spread = sapply(res, function(x) x[[1]]["spread"]),
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
  empty <- sapply(res, function(x) length(x) == 0)
  if(any(empty)) {
    res <- res[-which(empty)]
  }
  nclusts <- length(res)
  # print(nclusts)
  results <- NULL
  try(results <- data.frame(id = runif(1),
                        clust = 1:length(res),
                        size = sapply(res, function(x) x[[1]]["size"]),
                        snr = sapply(res, function(x) x[[1]]["snr"]),
                        rho = sapply(res, function(x) x[[1]]["rho"]),
                        threshold = sapply(res, function(x) x[[1]]["BHlevel"]),
                        spread = sapply(res, function(x) x[[1]]["spread"]),
                        slack = sapply(res, function(x) x[[1]]["slack"]),
                        true = sapply(res, function(x) x[[3]]["true"]),
                        cover = sapply(res, function(x) x[[2]][2] >= 0.05)))
  if(is.null(results)) {
    return(NULL)
  }

  return(results)
}
cover <- lapply(results, compCover)
cover <- data.table::rbindlist(cover)
cover <- data.frame(cover)

sumcover <- group_by(cover, id, rho, snr, threshold, spread, slack) %>%
  # summarize(cover = weighted.mean(cover, size))
  summarize(cover = mean(cover))
sumcover <- group_by(sumcover, rho, threshold, snr, spread, slack) %>%
  summarize(coversd = sd(cover) / sqrt(length(cover)),
            cover = mean(cover)) %>% as.data.frame()

ciquant <- qnorm(1 - 0.05/nrow(sumcover))
ggplot(sumcover, aes(x = snr, y = cover, col = factor(slack), linetype = factor(slack))) +
  geom_line() +
  geom_segment(linetype = 2, aes(xend = snr, y = cover - ciquant*coversd, yend = pmin(cover + ciquant*coversd, 1))) +
  facet_grid(rho ~ threshold, labeller = "label_both") +
  geom_hline(yintercept = 0.95) +
  theme_bw()

# where are we messing up ?
# ggplot(cover, aes(x = size, y = as.numeric(cover), col = factor(rho))) + geom_smooth() +
#   facet_grid(. ~ threshold, labeller = "label_both", scales = "free") + theme_bw()
#
# ggplot(cover, aes(x = true, y = as.numeric(cover), col = factor(rho))) + geom_smooth() +
#   facet_grid(. ~ threshold, labeller = "label_both", scales = "free") + theme_bw()

# Profile error -----------
# profileEstimates <- function(res) {
#   empty <- sapply(res, function(x) length(x) == 0)
#   if(any(empty)) {
#     res <- res[-which(empty)]
#   }
#   iterRes <- NULL
#   try(iterRes <- data.frame(id = runif(1),
#                             clust = 1:length(res),
#                             size = sapply(res, function(x) x[[1]]["size"]),
#                             snr = sapply(res, function(x) x[[1]]["snr"]),
#                             rho = sapply(res, function(x) x[[1]]["rho"]),
#                             threshold = sapply(res, function(x) x[[1]]["BHlevel"]),
#                             spread = sapply(res, function(x) x[[1]]["spread"]),
#                             slack = sapply(res, function(x) x[[1]]["slack"]),
#                             true = sapply(res, function(x) mean(x[[4]]$signal[x[[5]]])),
#                             prof = sapply(res, function(x) mean(x[[4]]$profile[x[[5]]]))))
#   return(iterRes)
# }
# profEsts <- lapply(results, profileEstimates)
# profEsts <- data.table::rbindlist(profEsts)
# profEsts <- data.frame(profEsts)
# profEsts$diff <- abs(profEsts$true - profEsts$prof)
# profEsts[order(profEsts$diff, decreasing = TRUE), ]
#
# res <- results[[sample.int(length(results), 1)]]
# res <- res[[sample.int(length(res), 1)]]
# plot(res[[4]]$profile, res[[4]]$signal, col = res[[5]] + 1,
#      main = paste("True Signal:", round(mean(res[[4]]$signal[res[[5]]]), 4)))
# abline(a = 0, b = 1)


# estimation error --------
library(reshape2)
compEstError <- function(res) {
  nclusts <- length(res)
  # print(nclusts)
  results <- NULL
  try(results <- data.frame(id = runif(1),
                            clust = 1:length(res),
                            size = sapply(res, function(x) x[[1]]["size"]),
                            snr = sapply(res, function(x) x[[1]]["snr"]),
                            rho = sapply(res, function(x) x[[1]]["rho"]),
                            threshold = sapply(res, function(x) x[[1]]["BHlevel"]),
                            spread = sapply(res, function(x) x[[1]]["spread"]),
                            slack = sapply(res, function(x) x[[1]]["slack"]),
                            conditional = sapply(res, function(x) (x[[3]]["conditional"] - x[[3]]["true"]) * sign(x[[3]]["naive"])),
                            naive = sapply(res, function(x) (x[[3]]["naive"] - x[[3]]["true"]) * sign(x[[3]]["naive"]))))
  if(is.null(results)) {
    return(NULL)
  }

  return(results)
}
estErr <- lapply(results, compEstError) %>%
  data.table::rbindlist() %>% data.frame()
sumEst <- melt(estErr, id = c("id", "clust", "size", "snr", "rho", "threshold", "spread", "slack"))
names(sumEst)[9:10] <- c("estimate", "error")
sumEst <- group_by(sumEst, id, snr, rho, threshold, spread, estimate, slack) %>%
  summarize(error = weighted.mean(error, size))

ggplot(subset(sumEst)) + geom_boxplot(aes(x = factor(snr), y = error, col = estimate)) +
  facet_grid(threshold ~ rho, labeller = "label_both") +
  theme_bw() + geom_hline(yintercept = 0)

ggplot(subset(sumEst)) + geom_boxplot(aes(x = factor(snr), y = abs(error), col = factor(estimate))) +
  facet_grid(threshold ~ rho, labeller = "label_both") +
  theme_bw() + geom_hline(yintercept = 0) + ylab("RMSE")







