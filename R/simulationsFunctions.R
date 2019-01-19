#' @export
findClusters <- function(coordinates, metric = "manhattan", distThreshold = 2) {
  columns <- which(names(coordinates) %in% c("i", "j", "k"))
  selected <- coordinates[coordinates$selected, ]
  graph <- as.matrix(dist(coordinates[, columns], method = metric))
  graph[graph > distThreshold] <- 0
  clusterNumber <- 1
  clusters <- list()
  while(nrow(selected) > 0) {
    cluster <- selected[1, columns]
    toadd <- which(graph[selected$row[1], ] != 0)
    toVerify <- toadd[which(coordinates$selected[toadd])]
    cluster <- rbind(cluster, coordinates[setdiff(toadd, toVerify), columns])
    selected <- selected[-1, ]
    while(length(toVerify) > 0) {
      srow <- which(selected$row == toVerify[1])
      toVerify <- toVerify[-1]
      cluster <- rbind(cluster, selected[srow, columns])
      toadd <- which(graph[selected$row[srow], ] != 0)
      newVerify <- toadd[which(toadd %in% setdiff(selected$row, toVerify))]
      toVerify <- c(toVerify, newVerify)
      selected <- selected[-srow, ]
      cluster <- rbind(cluster, coordinates[setdiff(toadd, toVerify), columns])
    }

    cluster <- unique(cluster)
    cluster$row <- as.numeric(rownames(cluster))
    cluster$selected <- coordinates$selected[cluster$row]
    clusters[[clusterNumber]] <- cluster
    clusterNumber <- clusterNumber + 1
  }

  return(clusters)
}

#' @export
generateArrayData3D <- function(dims, sqrtCov, targetSnr, spread = 1) {
  # parameters + setup
  I <- dims[1]
  J <- dims[2]
  K <- dims[3]

  # Generating Signal ------------
  coordinates <- expand.grid(i = 1:I, j = 1:J, k = 1:K)
  nnodes <- I * J * K
  coordinates$row <- 1:nrow(coordinates)

  # mu <- sapply(c(I, J, K), function(x) rnorm(x))
  # mu <- apply(coordinates, 1, function(x) {
  #   sum(mu[[1]][1:x[1]]) + sum(mu[[2]][x[2]:J]) + sum(mu[[3]][1:x[3]])
  # })
  # mu <- mu - mean(mu)
  location <- sapply(c(I, J, K), function(x) sample.int(x, 1))
  s <- matrix(0.3, nrow = 3, ncol = 3)
  diag(s) <- 1
  s <- s * spread
  mu <- mvtnorm::dmvnorm(coordinates[, 1:3], mean = location, sigma = s)
  mu <- mu * targetSnr / max(mu)
  coordinates$signal <- mu

  # Generating noise + data -----------------
  noise <- rnorm(nnodes)
  noise <- sqrtCov %*% noise
  coordinates$noise <- noise
  coordinates$observed <- coordinates$signal + coordinates$noise

  return(coordinates)
}

#' @export
residualData3D <- function(dat, grp_size, targetSnr ,spread = 1) {

  dims = dim(dat)

  # Generating signal
  coordinates <- expand.grid(i=1:dims[1],j=1:dims[2],k=1:dims[3])
  nnodes <- prod(dims[1:3])
  coordinates$row = 1:nrow(coordinates)

  location <- sapply(dims[1:3], function(x) sample.int(x,1))
  s <- matrix(0.3, nrow = 3, ncol = 3)
  diag(s) <- 1
  s <- s * spread
  mu <- mvtnorm::dmvnorm(coordinates[, 1:3], mean = location,
                         sigma = s)
  mu <- mu * targetSnr/max(mu)
  coordinates$signal <- round(mu, 3)

  # Generating noise + data -----------------
  pop_n = dim(dat)[4]
  samp = sample(pop_n, grp_size*2, replace = FALSE)
  noise <- as.numeric(apply(dat[,,,samp[1:grp_size]],c(1,2,3),mean) -
                        apply(dat[,,,samp[grp_size+(1:grp_size)]],c(1,2,3),mean))

  scale_noise = sd(noise)
  coordinates$noise <- noise/scale_noise
  coordinates$scale_coef = scale_noise
  coordinates$observed <- coordinates$signal + coordinates$noise

  return(coordinates)
}

