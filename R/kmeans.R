#' @title K-Means Clustering
#'
#' @description Perform K-Means clustering on a data matrix.
#'
#' @param data a set of observations, presented as a matrix-like object where
#' every row is a new observation.
#' @param centers either the number of clusters or a set of initial cluster
#' centers. If a number, the centers are chosen according to the
#' \code{initialization} parameter.
#' @param max_iterations the maximum number of iterations allowed.
#' @param initialization the initialization method to be used. This should be
#' one of \code{"random"} or \code{"kmeans++"}. The latter is the default.
#' @param ... additional arguments passed to [proxy::dist()].
#'
#' @details The data given by \code{data} is clustered by the \eqn{k}-means
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
#'    \item Computation of the distance from each data point to each center.
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
#' The \code{initialization} methods provided by this function are:
#'
#' \describe{
#'  \item{\code{random}}{A set of \code{centers} observations is chosen at
#'  random from the data as the initial centers.}
#'  \item{\code{kmeans++}}{The \code{centers} observations are chosen using the
#'  \strong{kmeans++} algorithm. This algorithm chooses the first center at
#'  random and then chooses the next center from the remaining observations with
#'  probability proportional to the square distance to the closest center. This
#'  process is repeated until \code{centers} centers are chosen.}
#' }
#'
#' @return A [stats::kmeans()] object.
#'
#' @examples
#' ### Voronoi tesselation
#' voronoi <- suppressMessages(suppressWarnings(require(deldir)))
#' cols <- c(
#'   "#00000019",
#'   "#DF536B19",
#'   "#61D04F19",
#'   "#2297E619",
#'   "#28E2E519",
#'   "#CD0BBC19",
#'   "#F5C71019",
#'   "#9E9E9E19"
#' )
#'
#' ### Helper function
#' test <- function(db, k) {
#'   print(cl <- clustlearn::kmeans(db, k, 100))
#'   plot(db, col = cl$cluster, asp = 1, pch = 20)
#'   points(cl$centers, col = seq_len(k), pch = 13, cex = 2, lwd = 2)
#'
#'   if (voronoi) {
#'     x <- c(min(db[, 1]), max(db[, 1]))
#'     dx <- c(x[1] - x[2], x[2] - x[1])
#'     y <- c(min(db[, 2]), max(db[, 2]))
#'     dy <- c(y[1] - y[2], y[2] - y[1])
#'     tesselation <- deldir(
#'       cl$centers[, 1],
#'       cl$centers[, 2],
#'       rw = c(x + dx, y + dy)
#'     )
#'     tiles <- tile.list(tesselation)
#'
#'     plot(
#'       tiles,
#'       asp = 1,
#'       add = TRUE,
#'       showpoints = FALSE,
#'       border = "#00000000",
#'       fillcol = cols
#'     )
#'   }
#' }
#'
#' ### Example 1
#' test(clustlearn::db1, 2)
#'
#' ### Example 2
#' test(clustlearn::db2, 2)
#'
#' ### Example 3
#' test(clustlearn::db3, 3)
#'
#' ### Example 4
#' test(clustlearn::db4, 3)
#'
#' ### Example 5
#' test(clustlearn::db5, 3)
#'
#' ### Example 6
#' test(clustlearn::db6, 3)
#'
#' @author Eduardo Ruiz Sabajanes, \email{eduardoruizsabajanes@@gmail.com}
#'
#' @importFrom proxy dist
#' @export
kmeans <- function(
  data,
  centers,
  max_iterations = 10,
  initialization = "kmeans++",
  ...
) {
  # Make sure max_iterations is a positive integer
  if (!is.numeric(max_iterations) || max_iterations < 1)
    stop("max_iterations must be an integer greater than 0")

  # Get centers
  if (missing(centers))
    stop("centers must be a matrix or a number")

  if (length(centers) == 1) {
    if (centers < 1)
      stop("centers must be a positive integer")
    if (centers > nrow(data))
      stop("centers must be less than or equal to the number of observations")

    # Figure out the initialization method
    initialization <- grep(
      tolower(initialization),
      c("random", "kmeans++"),
      fixed = TRUE
    )

    if (length(initialization) != 1)
      stop("initialization must be one of 'random' or 'kmeans++'")

    # Initialize centers ...
    centers <- switch(
      initialization,

      # ... randomly
      random_init(as.matrix(data), centers, ...),

      # ... using the kmeans++ algorithm
      kmeanspp_init(as.matrix(data), centers, ...)
    )
  }

  # Update centers while they don't change
  iter <- 0
  for (i in seq_len(max_iterations)) {
    iter <- i
    old_centers <- as.matrix(centers)

    # Compute distances between points and centers
    distances <- proxy::dist(old_centers, data, ...)

    # Find which center is closest to each point
    nearest_centers <- apply(distances, 2, which.min)

    # Compute the new centers as the average of it's closest points
    new_centers <- matrix(
      sapply(
        seq_len(nrow(old_centers)),
        function(n) {
          temp <- data[nearest_centers == n, , drop = FALSE]
          if (nrow(temp) > 0)
            apply(temp, 2, mean)
          else
            old_centers[n, ]
        }
      ),
      nrow = ncol(old_centers)
    )
    centers <- t(new_centers)

    # If centers aren't updated
    if (all(centers == old_centers))
      break
  }

  # Compute distances between points and centers
  distances <- proxy::dist(centers, data, ...)

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
  tot.withinss <- sum(withinss)

  # The between-cluster sum of squares
  betweenss <- totss - tot.withinss

  # Find the size of each cluster
  tmp <- factor(nearest_centers, levels = seq_len(nrow(centers)))
  size <- as.integer(table(tmp))

  structure(
    list(
      cluster = nearest_centers,
      centers = centers,
      totss = totss,
      withinss = withinss,
      tot.withinss = tot.withinss,
      betweenss = betweenss,
      size = size,
      iter = iter,
      ifault = 0
    ),
    class = "kmeans"
  )
}

random_init <- function(data, k, ...) {
  smp <- sample(nrow(data), size = k, replace = FALSE)
  data[smp, , drop = FALSE]
}

kmeanspp_init <- function(data, k, ...) {
  centers <- matrix(0, nrow = k, ncol = ncol(data))
  probs <- rep(1 / nrow(data), nrow(data))
  for (i in seq_len(k)) {
    # Choose a center with probability proportional to its square distance to
    # the closest center
    smp <- sample(nrow(data), size = 1, replace = FALSE, prob = probs)
    centers[i, ] <- data[smp, ]

    # Update the probabilities
    distances <- proxy::dist(centers[seq_len(i), , drop = FALSE], data, ...) ^ 2
    probs <- apply(distances, 2, min)
    probs <- probs / sum(probs)

    # Replace NAs and NaNs with the remaining probability
    tmp <- sum(probs[is.finite(probs)])
    probs[!is.finite(probs)] <- (1 - tmp) / sum(!is.finite(probs))
  }
  centers
}
