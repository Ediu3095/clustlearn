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
#' @param ... additional arguments passed to [proxy::dist()].
#'
#' @details The data given by \code{data} are clustered by the \eqn{k}-means
#' method, which aims to partition the points into \eqn{k} groups such that the
#' sum of squares from points to the assigned cluster centers is minimized. At
#' the minimum, all cluster centers are at the mean of their Voronoi sets (the
#' set of data points which are nearest to the cluster center).
#'
#' The \eqn{k}-means method follows a 2 to \eqn{n} step process:
#'
#' \enumerate{
#'  \item The first step can be subdivided into 3 steps: \enumerate{
#'    \item Selection of the number \eqn{k} of clusters, into which the data is
#'    going to be grouped and of which the centers will be the representatives.
#'    This is determined through the use of the \code{centers} parameter.
#'    \item Computation of the euclidean distance to each of the centers.
#'    \item Assignment of each observation to a cluster. The observation is
#'    assigned to the cluster represented by the nearest center.
#'  }
#'  \item The next steps are just like the first but for the first sub-step:
#'  \enumerate{
#'    \item Computation of the new centers. The center of each cluster is
#'    computed as the mean of the observations assigned to said cluster.
#'  }
#' }
#'
#' The algorithm stops once the centers in step \eqn{n+1} are the same as the
#' ones in step \eqn{n}. However, this convergence does not always take place.
#' For this reason, the algorithm also stops once a maximum number of iterations
#' \code{max_iterations} is reached.
#'
#' @author Eduardo Ruiz Sabajanes, \email{eduardo.ruizs@@edu.uah.es}
#'
#' @return A [stats::kmeans()] object.
#'
#' @examples
#' # a 2-dimensional example
#' x <- rbind(matrix(rnorm(100, sd = 0.3), ncol = 2),
#'            matrix(rnorm(100, mean = 1, sd = 0.3), ncol = 2))
#' colnames(x) <- c("x", "y")
#' (cl <- kmeans(x, 2))
#' plot(x, col = cl$cluster, pch = 20)
#' points(cl$centers, col = 1:2, pch = 10, cex = 2)
#'
#' # sum of squares
#' ss <- function(x) sum(scale(x, scale = FALSE)^2)
#'
#' ## cluster centers "fitted" to each obs.:
#' fitted.x <- fitted(cl);  head(fitted.x)
#' resid.x <- x - fitted(cl)
#'
#' @importFrom stats runif
#' @export
kmeans <- function(data, centers, max_iterations = 10, ...) {
  # Get centers
  if (missing(centers))
    stop('centers must be a matrix or a number')
  if (length(centers) == 1) {
    smp <- sample(nrow(data), size = centers, replace = FALSE)
    centers <- data[smp, , drop = FALSE]
  }

  # Compute distances between points and centers
  distances <- proxy::dist(data, centers, ...)
  iter <- 0

  # Update centers while they don't change
  for (i in seq_len(max_iterations)) {
    iter <- i
    old_centers <- as.matrix(centers)

    # Compute distances between points and centers
    distances <- proxy::dist(data, old_centers, ...)

    # Find which center is closest to each point
    nearest_centers <- apply(distances, 2, which.min)

    # Compute the new centers as the average of it's closest points
    new_centers <- sapply(
      seq_len(nrow(old_centers)),
      function(n) {
        temp <- data[nearest_centers == n, , drop = FALSE]
        if (nrow(temp) > 0)
          apply(temp, 2, mean)
        else
          old_centers[n, ]
      }
    )
    centers <- t(new_centers)

    # If centers aren't updated
    if (all(centers == old_centers))
      break
  }

  # Find which center is closest to each point
  nearest_centers <- apply(distances, 2, which.min)

  # Name the centroids
  row.names(centers) <- seq_len(nrow(centers))

  # Total sum of squares
  center <- apply(data, 2, mean)
  totss <- sum(apply(data, 1, function(x) x - center)^2)

  # Total within-cluster sum of squares
  withinss <- sapply(
    seq_len(nrow(centers)),
    function(cluster) {
      ccenter <- centers[cluster, ]
      cdata <- data[nearest_centers == cluster, , drop = FALSE]
      sum(apply(cdata, 1, function(x) x - ccenter)^2)
    }
  )

  # Total within-cluster sum of squares
  tot_withinss <- sum(withinss)

  # The between-cluster sum of squares
  betweenss <- totss - tot_withinss

  structure(
    list(
      cluster = nearest_centers,
      centers = centers,
      totss = totss,
      withinss = withinss,
      tot.withinss = tot_withinss,
      betweenss = betweenss,
      size = as.integer(table(nearest_centers)),
      iter = iter,
      ifault = 0
    ),
    class = "kmeans"
  )
}
