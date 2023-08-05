#' @title Density Based Spatial Clustering of Applications with Noise (DBSCAN)
#'
#' @description Cluster analysis to find clusters of arbitrary shape and is
#' robust to outliers.
#'
#' @param data a set of observations, presented as a matrix-like object where
#' every row is a new observation.
#' @param eps how close two observations have to be to be considered neighbors.
#' @param min_pts the minimum amount of neighbors for a region to be considered
#' dense.
#' @param ... additional arguments passed to [proxy::dist()].
#'
#' @details
#'
#' @author Eduardo Ruiz Sabajanes, \email{eduardo.ruizs@@edu.uah.es}
#'
#' @return A [dbscan()] object.
#'
#' @examples
#' # a 2-dimensional example
#' x <- rbind(matrix(rnorm(100, sd = 0.3), ncol = 2) + c(0, 0),
#'            matrix(rnorm(100, sd = 0.3), ncol = 2) + c(1, 1))
#' colnames(x) <- c("x", "y")
#' (cl <- dbscan(x, 0.1, 5))
#' plot(x, col = cl + 1, pch = 20)
#'
#' @importFrom proxy dist
#' @export
dbscan <- function(
  data,
  eps = .5,
  min_pts = 5,
  ...
) {
  # Precompute neighbors
  distances <- proxy::dist(data, ...)
  neighbors <- distances <= eps

  # Initialize clusters
  cluster_id <- new.env()
  cluster_id$current <- 1
  cluster_id$of <- rep(-1, nrow(data))

  # Each loop finds a new cluster around a core point
  for (idx in seq_len(nrow(data))) {
    if (cluster_id$of[idx] != -1)
      next
    if (expand_cluster(neighbors, cluster_id, idx, min_pts))
      cluster_id$current <- cluster_id$current + 1
  }

  # Return a dbscan object
  structure(
    list(
      cluster = cluster_id$of,
      eps = eps,
      min_pts = min_pts
    ),
    class = "dbscan"
  )
}

expand_cluster <- function(
  neighbors,
  cluster_id,
  point,
  min_pts
) {
  # Get the point's neighbors (including itself)
  seeds <- region_query(neighbors, point)

  if (length(seeds) < min_pts) {
    # If it is not a core point, it is noise
    cluster_id$of[point] <- 0
    FALSE
  } else {
    # Otherwise, we can expand the cluster
    cluster_id$of[seeds] <- cluster_id$current
    frontier <- setdiff(seeds, point)

    # Loop until there are no more neighbors in the frontier
    while (length(frontier) > 0) {
      current_point <- frontier[1]

      # Get current_point's neighbors
      result <- region_query(neighbors, current_point)

      # If current_point is a core point, expand the cluster
      if (length(result) >= min_pts) {
        # Add non visited neighbors to the frontier
        not_visited <- cluster_id$of[result] == -1
        frontier <- c(frontier, result[not_visited])

        # Add not clustered neighbors to the cluster
        noise <- cluster_id$of[result] == 0
        cluster_id$of[result][not_visited | noise] <- cluster_id$current
      }

      # Remove current_point from the frontier
      frontier <- frontier[-1]
    }
    TRUE
  }
}

region_query <- function(
  neighbors,
  idx
) {
  # Return the indices of the neighbors of idx
  which(neighbors[idx, ])
}
