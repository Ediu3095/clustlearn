#' @title Divisive Hierarchical Clustering
#'
#' @description Hierarchical divisive cluster analysis on a set of observations
#'
#' @param data a set of observations, presented as a matrix-like object where
#' every row is a new observation.
#' @param eps how close two observations have to be to be considered neighbors.
#' @param min_pts the minimum amount of neighbors for a region to be considered
#' dense.
#' @param method the distance definition to be used. Check [distance()] for the
#' available methods.
#' @param p the power of the Minkowski distance to be used. Only used if
#' \code{method = "minkowski"}.
#'
#' @details
#'
#' @author Eduardo Ruiz Sabajanes, \email{eduardo.ruizs@@edu.uah.es}
#'
#' @return A [kmeans()] object.
#'
#' @examples
#' # a 2-dimensional example
#' x <- rbind(matrix(rnorm(100, sd = 0.3), ncol = 2),
#'            matrix(rnorm(100, mean = 1, sd = 0.3), ncol = 2))
#' colnames(x) <- c("x", "y")
#' (cl <- dbscan(x, 0.5, 5))
#' plot(x, col = cl, pch = 20)
#'
#' @importFrom stats dist
#' @export
dbscan <- function(
  data,
  eps = .5,
  min_pts = 5,
  method = "euclidean",
  p = 2
) {
  distances <- distance(data, method = method, p = p, simplify = FALSE)
  thresholds <- distances <= eps

  # A point is dense if it has at least min_pts neighbors within eps
  is_dense <- function(i) sum(distances[i, ] <= eps) > min_pts

  # Identify dense regions
  dense_regions <- sapply(seq_len(nrow(data)), is_dense)

  # Initialize clusters
  clusters <- rep(-1, nrow(data))

  # Elegible points are those that are dense and have not been clustered yet
  elegible_points <- dense_regions & (clusters == -1)

  # Each loop finds a new cluster around a core point
  while (any(elegible_points)) {
    # Pick a random core point
    core_point <- sample(which(elegible_points), 1)

    # Assign a new cluster to the core point
    clusters[core_point] <- max(clusters) + 1

    # Find all core points within eps of the selected core point or any of its
    # neighbors
    neighbors <- NULL
    new_neighbors <- core_point
    while (length(new_neighbors) > 0) {
      neighbors <- c(neighbors, new_neighbors)
      new_neighbors <- thresholds[new_neighbors, , drop = FALSE]
      new_neighbors <- apply(new_neighbors, 2, any)
      new_neighbors <- new_neighbors & dense_regions
      new_neighbors <- which(new_neighbors)
      new_neighbors <- setdiff(new_neighbors, neighbors)
    }
    # Find all non-core points within eps of the selected core point or any of
    # its neighbors
    new_neighbors <- thresholds[neighbors, , drop = FALSE]
    new_neighbors <- apply(new_neighbors, 2, any)
    new_neighbors <- which(new_neighbors)
    neighbors <- new_neighbors

    # Assign the same cluster to all neighbors
    clusters[neighbors] <- clusters[core_point]

    # Update elegible points
    elegible_points <- dense_regions & (clusters == -1)
  }

  clusters <- clusters + 1

  clusters
}
