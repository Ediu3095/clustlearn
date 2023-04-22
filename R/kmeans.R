#' @title K-Means Clustering
#'
#' @description Perform K-Means clustering on a data matrix.
#'
#' @param data a set of observations, presented as a matrix-like object where
#' every row is a new observation.
#' @param centers either the number of clusters or a set of initial cluster
#' centers. If a number, a random set of rows in x is chosen as the initial
#' centers.
#' @param max_iterations the maximum number of iterations allowed.
#' @param method the distance definition to be used. Check [distance()] for the
#' available methods.
#' @param p the exponent of the Minkowski distance.
#'
#' @details
#'
#' @author Eduardo Ruiz Sabajanes, \email{eduardo.ruizs@@edu.uah.es}
#'
#' @return A [kmeans()] object.
#'
#' @examples
#'
#' @importFrom stats runif
#' @export
kmeans <- function(data, centers, max_iterations = 10, method = "euclidean", p = 3) {
  # Get centers
  if (missing(centers))
    stop('centers must be a matrix or a number')
  if (length(centers) == 1)
    centers <- apply(
      data,
      2,
      function(col) {
        runif(centers, min = min(col), max = max(col))
      })

  # Update centers while they don't change
  for (i in seq_len(max_iterations)) {
    iter <- i
    old_centers <- centers

    # Compute distances between points and centers
    distances <- distance(data, old_centers)

    # Find which center is closest to each point
    nearest_centers <- apply(distances, 2, which.min)

    # Compute the new centers as the average of it's closest points
    new_centers <- sapply(
      seq_len(nrow(old_centers)),
      function(n) {
        temp <- data[nearest_centers == n, ]
        if (nrow(temp) > 0)
          apply(temp, 2, mean)
        else
          old_centers[n, ]
      })
    centers <- t(new_centers)

    # If centers aren't updated
    if (all(centers == old_centers))
      break
  }

  # Find which center is closest to each point
  nearest_centers <- apply(distances, 2, which.min)

  # Name the centroids
  row.names(centers) <- seq_len(nrow(centers))

  structure(
    list(
      cluster = nearest_centers,
      centers = centers,
      totss = NULL,
      withinss = NULL,
      tot.withinss = NULL,
      betweenss = NULL,
      size = as.integer(table(nearest_centers)),
      iter = iter,
      ifault = 0
    ),
    class = 'kmeans'
  )
}
