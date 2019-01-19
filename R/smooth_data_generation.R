#' @import mvtnorm
#' @importFrom  expm sqrtm
#' @import magrittr
generate_AR_grid_noise <- function(n, dims, noise_sd, rho = 0.0,
                                   dist_method = "manhattan") {
  if(any(dims < 1)) {
    stop("All dimensions must be >= 1")
  }

  coordinates <- expand.grid(lapply(dims, function(x) 1:x))
  distances <- as.matrix(dist(coordinates, method = dist_method))
  noise_covariance <- rho^distances * noise_sd^2
  noise <- rmvnorm(n, sigma = noise_covariance)
  return(list(noise = noise,
              coordinates = coordinates,
              covariance = noise_covariance,
              distances = distances))
}

generate_smooth_signal <- function(distances, mu_sd, bandwidth_sd) {
  n_coordinates <- nrow(distances)
  mu <- rnorm(n_coordinates, sd = mu_sd)
  smooth_mu <- numeric(n_coordinates)
  for(i in 1:n_coordinates) {
    weights <- dnorm(distances[i, ], sd = bandwidth_sd)
    smooth_mu[i] <- weighted.mean(mu, weights)
  }
  return(smooth_mu)
}

#' @export
generate_smooth_coordinate_data <- function(n, dims, rho = 0.0,
                                            noise_sd = 1,
                                            mu_sd = 1,
                                            bandwidth_sd = 3,
                                            dist_method = "manhattan") {
  noise <- generate_AR_grid_noise(n, dims, noise_sd, rho, dist_method)
  distances <- noise$distances
  covariance <- noise$covariance
  coordinates <- noise$coordinates
  noise <- noise$noise
  mu <- replicate(n, generate_smooth_signal(distances, mu_sd, bandwidth_sd)) %>% t()
  return(list(y = noise + mu, mu = mu, noise = noise,
              covariance = covariance,
              coordinates = coordinates,
              distances = distances))
}

#' find clusters of adjascent selected coordinates
#'
#' @param selected_coordinates an integer vector, the coordinates of which
#' denotes which coordinates have crossed the threshold (have been selected)
#' such that will be obtained, for example from calling
#' \code{which(x > threshold)}.
#'
#' @export
find_clusters <- function(distances, selected_coordinates,
                          magnitude_threshold,
                          dist_threshold = 2, dist_method = "manhattan") {
  n_coordinates <- nrow(distances)
  cluster_list <- vector(n_coordinates, mode = "list")
  list_index <- 1
  needs_assignment <- rep(TRUE, length(selected_coordinates))
  while(sum(needs_assignment) > 0) {
    to_check <- selected_coordinates[needs_assignment][1]
    cluster <- numeric(sum(needs_assignment))
    cluster_index <- 1
    needs_assignment[needs_assignment][1] <- FALSE
    while(length(to_check) > 0) {
      new_candidate <- to_check[1]
      to_check <- to_check[-1]
      cluster[cluster_index] <- new_candidate
      cluster_index <- cluster_index + 1
      neighbors <- distances[new_candidate, selected_coordinates[needs_assignment]] <= dist_threshold
      to_check <- union(selected_coordinates[needs_assignment][neighbors], to_check)
      if(any(neighbors)) {
        needs_assignment[needs_assignment][neighbors] <- FALSE
      }
      #print(to_check)
    }
    cluster_list[[list_index]] <- cluster[1:(cluster_index - 1)]
    list_index <- list_index + 1
  }

  return(cluster_list[1:(list_index - 1)])
}


