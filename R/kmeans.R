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
#' @param details a Boolean determining whether intermediate logs explaining how
#' the algorithm works should be printed or not.
#' @param waiting a Boolean determining whether the intermediate logs should be
#' printed in chunks waiting for user input before printing the next or not.
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
#' ### Example 7 (with explanations, no plots)
#' cl <- clustlearn::kmeans(
#'   clustlearn::db5[1:20, ],
#'   3,
#'   details = TRUE,
#'   waiting = FALSE
#' )
#'
#' @author Eduardo Ruiz Sabajanes, \email{eduardo.ruizs@@edu.uah.es}
#'
#' @importFrom proxy dist
#' @export
kmeans <- function(
  data,
  centers,
  max_iterations = 10,
  initialization = "kmeans++",
  details = FALSE,
  waiting = TRUE,
  ...
) {
  # Make sure max_iterations is a positive integer
  if (!is.numeric(max_iterations) || max_iterations < 1)
    stop("max_iterations must be an integer greater than 0")

  # Get centers
  if (missing(centers))
    stop("centers must be a matrix or a number")

  init <- 3
  if (length(centers) == 1) {
    if (centers < 1)
      stop("centers must be a positive integer")
    if (centers > nrow(data))
      stop("centers must be less than or equal to the number of observations")

    # Figure out the initialization method
    init <- grep(
      tolower(initialization),
      c("random", "kmeans++"),
      fixed = TRUE
    )

    if (length(init) != 1)
      stop("initialization must be one of 'random' or 'kmeans++'")
  }

  if (details) {
    hline()
    console.log("EXPLANATION:")
    console.log("")
    console.log("The K-Means algorithm aims to partition a dataset into k groups such that the within-cluster sum-of-squares is minimized. At the minimum, all cluster centers are at the mean of their Voronoi sets (the set of data points which are nearest to the cluster center).")
    console.log("")
    console.log("The K-Means method follows a 2 to n step process:")
    console.log("")
    console.log("    1. The first step can be subdivided into 3 steps:")
    console.log("        1. Selection of the number k of clusters, into which the data is going to be grouped and of which the centroids will be the representatives.")
    console.log("        2. Computation of the distance from each observation to each centroid.")
    console.log("        3. Assignment of each observation to a cluster. Observations are assigned to the cluster represented by the nearest centroid.")
    console.log("    2. The next steps are just like the first but for the first sub-step we do the following:")
    console.log("        1. Computation of the new centroids The centroid of each cluster is computed as the mean of the observations assigned to said cluster.")
    console.log("")
    console.log("The algorithm stops once the centers in step n+1 are the same as the ones in step n. However, this convergence does not always take place. For this reason, the algorithm also stops once a maximum number of iterations is reached.")
    console.log("")

    if (waiting) {
      invisible(readline(prompt = "Press [enter] to continue"))
      console.log("")
    }

    hline()
    console.log("STEP: 1")
    console.log("")
    console.log("If they are not, k centroids have to be initialized...")
    console.log("")
  }

  # Initialize centers ...
  centers <- switch(
    init,
    # ... randomly
    random_init(as.matrix(data), centers, details, waiting, ...),
    # ... using the kmeans++ algorithm
    kmeanspp_init(as.matrix(data), centers, details, waiting, ...),
    # ... they are already initialized
    centers
  )

  if (details) {
    console.log("With this, the k initial centroids are the following:")
    cat("Centroids:\n")
    print(centers)
    console.log("")

    if (waiting) {
      invisible(readline(prompt = "Press [enter] to continue"))
      console.log("")
    }
  }

  # Update centers while they don't change
  iter <- 0
  for (i in seq_len(max_iterations)) {
    iter <- i
    old_centers <- as.matrix(centers)

    # Compute distances between points and centers
    distances <- proxy::dist(old_centers, data, ...)

    if (details) {
      console.log("With these centroids, the matrix of pairwise distances between the observations and the centroids is computed:")
      cat("Distances:\n")
      print(round(distances, 3))
      console.log("")

      if (waiting) {
        invisible(readline(prompt = "Press [enter] to continue"))
        console.log("")
      }
    }

    # Find which center is closest to each point
    nearest_centers <- apply(distances, 2, which.min)

    if (details) {
      console.log("From these distances, the centroid closest to each observation is computed. In this way, we make the following observation-cluster assignments:")
      cat("Cluster assignments:\n")
      print(nearest_centers)
      console.log("")

      if (waiting) {
        invisible(readline(prompt = "Press [enter] to continue"))
        console.log("")
      }
    }

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

    if (details) {
      hline()
      console.log(paste("STEP:", iter + 1))
      console.log("")
      console.log("With the previous cluster assignments, we compute the new centroids as the mean of the observations assigned to the corresponding cluster:")
      cat("Centroids:\n")
      print(centers)
      console.log("")

      if (waiting) {
        invisible(readline(prompt = "Press [enter] to continue"))
        console.log("")
      }
    }

    # If centers aren't updated
    if (all(centers == old_centers))
      break
  }

  # Compute distances between points and centers
  distances <- proxy::dist(centers, data, ...)

  if (details) {
    if (iter == max_iterations)
      console.log("Since we have reached the maximum amount of iterations, these are the last centroids we are going to compute.")
    else
      console.log("Since the centroids have not changed with regards to last step, these are the last centroids we are going to compute.")
    console.log("With these centroids, the matrix of pairwise distances between the observations and the centroids is computed one last time:")
    cat("Distances:\n")
    print(round(distances, 3))
    console.log("")

    if (waiting) {
      invisible(readline(prompt = "Press [enter] to continue"))
      console.log("")
    }
  }

  # Find which center is closest to each point
  nearest_centers <- apply(distances, 2, which.min)
  row.names(centers) <- seq_len(nrow(centers))

  if (details) {
    console.log("From these distances, observations are assigned the cluster of whichever centroid they are closest to:")
    cat("Cluster assignments:\n")
    print(nearest_centers)
    console.log("")

    if (waiting) {
      invisible(readline(prompt = "Press [enter] to continue"))
      console.log("")
    }

    hline()
    console.log("")
  }

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

random_init <- function(data, k, details, waiting, ...) {
  if (details) {
    console.log("In this case, the k centroids are chosen randomly from the observations...")
    console.log("")
  }

  smp <- sample(nrow(data), size = k, replace = FALSE)
  data[smp, , drop = FALSE]
}

kmeanspp_init <- function(data, k, details, waiting, ...) {
  if (details) {
    console.log("In this case, the k centroids are chosen using the kmeans++ algorithm...")
    console.log("")
  }

  centers <- matrix(0, nrow = k, ncol = ncol(data))
  probs <- rep(1 / nrow(data), nrow(data))
  for (i in seq_len(k)) {
    # Choose a center with probability proportional to its square distance to
    # the closest center
    smp <- sample(nrow(data), size = 1, replace = FALSE, prob = probs)
    centers[i, ] <- data[smp, ]

    if (details) {
      console.log(paste0("--- kmeans++ step #", i, " ---\n"))
      console.log("A centroid is chosen according to the following probabilities:")
      cat("Probs:\n")
      print(round(probs, 3))
      console.log("The chosen centroid is:")
      cat(paste0("Observation #", smp, ":\n"))
      print(centers[i, ])
      cat("With probability:\n")
      print(probs[[smp]])
      if (i < k)
        console.log("With this new centroid the probabilities are updated...")
      else
        console.log("With this new centroid we already have k centroids...")
      console.log("")

      if (waiting) {
        invisible(readline(prompt = "Press [enter] to continue"))
        console.log("")
      }
    }

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
