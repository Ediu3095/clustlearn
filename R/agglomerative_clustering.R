#' @title Agglomerative Hierarchical Clustering
#'
#' @description Hierarchical agglomerative cluster analysis on a set of
#' observations
#'
#' @param data a set of observations, presented as a matrix-like object where
#' every row is a new observation.
#' @param proximity the proximity definition to be used. This should be one
#' of \code{"MIN"} (single linkage), \code{"MAX"} (complete linkage),
#' \code{"AVG"} (average linkage).
#' @param ... additional arguments passed to [proxy::dist()].
#'
#' @details This function performs a hierarchical cluster analysis for the
#' \eqn{n} objects being clustered. The definition of a set of clusters using
#' this method follows a \eqn{n} step process, which repeats until a single
#' cluster remains:
#'
#' \enumerate{
#'  \item Initially, each object is assigned to its own cluster. The matrix
#'  of distances between clusters is computed. This is done according to the
#'  specified \code{method}.
#'  \item The two clusters with closest proximity be joined together and the
#'  proximity matrix will be updated. This is done according to the specified
#'  \code{proximity}. This step is repeated until a single cluster remains.
#' }
#'
#' The definitions of \code{proximity} considered by this function are:
#'
#' \describe{
#'  \item{\code{MIN}}{\eqn{\min\left\{d(x,y):x\in A,y\in B\right\}}. Defines the
#'  proximity between two clusters as the distance between the closest objects
#'  among the two clusters. It produces clusters where each object is closest to
#'  at least one other object in the same cluster. It is also known as
#'  \strong{SLINK} or \strong{single-link}.}
#'  \item{\code{MAX}}{\eqn{\max\left\{d(x,y):x\in A,y\in B\right\}}. Defines the
#'  proximity between two clusters as the distance between the furthest objects
#'  among the two clusters. It is also known as \strong{CLINK} or
#'  \strong{complete-link}.}
#'  \item{\code{AVG}}{\eqn{\frac{1}{\left|A\right|\cdot\left|B\right|}
#'  \sum_{x\in A}\sum_{y\in B} d(x,y)}. Defines the proximity between two
#'  clusters as the average distance between every pair of objects, one from
#'  each cluster. It is also known as \strong{UPGMA} or \strong{average-link}.}
#' }
#'
#' @author Eduardo Ruiz Sabajanes, \email{eduardo.ruizs@@edu.uah.es}
#'
#' @return An [hclust()] object which describes the tree produced by the
#' clustering process.
#'
#' @examples
#' ### Example 1: Violent crime rates by US state
#'
#' ac <- agglomerative_clustering(USArrests, proximity = "AVG")
#' plot(ac, hang = -1)
#'
#' ### Example 2: Iris flower dataset
#'
#' ac <- agglomerative_clustering(iris[, 1:4], proximity = "AVG")
#'
#' # Dendrogram
#' plot(ac, hang = -1, label = iris$Species)
#' rect.hclust(ac, k = 3, border = 1:3)
#'
#' # Scatterplot
#' plot(iris, col = cutree(ac, k = 3))
#'
#' @importFrom proxy dist
#' @export
agglomerative_clustering <- function(data, proximity = "MIN", ...) {
  # Function needed to calculate the avg distance between two clusters
  avg <- function(m1, m2) function(d1, d2) (d1 * m1 + d2 * m2) / (m1 + m2)

  # Prepare the data structure which will hold the answer
  ans <- structure(
    list(
      merge = numeric(0),
      height = NULL,
      order = NULL,
      labels = rownames(data),
      method = proximity,
      call = NULL,
      dist.method = "euclidean"
    ),
    class = "hclust"
  )

  # Compute the distances between each point
  d <- as.matrix(proxy::dist(data, ...))
  d <- mapply(
    "[<-",
    data.frame(d),
    seq_len(nrow(data)),
    sample(Inf, nrow(data), TRUE),
    USE.NAMES = FALSE
  )
  method <- attr(d, "method")
  ans$dist.method <- if (is.null(method)) "Euclidean" else method

  # Create a list with the initial clusters
  c <- lapply(
    seq_len(nrow(data)),
    function(data) {
      structure(
        data,
        label = -data,
        members = 1
      )
    }
  )

  for (i in seq_len(length(c) - 1)) {
    # Look for the minimum distance between two clusters
    md <- which.min(d) - 1
    md <- sort(c(md %% nrow(d), md %/% nrow(d)) + 1)

    # Join the clusters into a new one
    c1 <- c[[md[1]]]
    m1 <- attr(c1, "members")
    c2 <- c[[md[2]]]
    m2 <- attr(c2, "members")
    c3 <- structure(
      list(c1, c2),
      label = i,
      members = m1 + m2
    )

    # Add the merged clusters to the answer
    ans$merge <- c(ans$merge, attr(c1, "label"), attr(c2, "label"))
    ans$height <- c(ans$height, d[md[1], md[2]])
    c[[md[1]]] <- c3
    c <- c[-md[2]]

    # Recompute the distances (proximity)
    d1 <- d[, md[1]]
    d2 <- d[, md[2]]
    d3 <- mapply(
      switch(
        proximity,
        MAX = max,
        MIN = min,
        AVG = avg(m1, m2)
      ),
      d1,
      d2
    )
    d3[md] <- Inf
    d[, md[1]] <- d3
    d[md[1], ] <- d3
    d <- d[-md[2], -md[2]]
  }

  # Compute the merge and order of the hclust
  ans$merge <- matrix(ans$merge, ncol = 2, byrow = TRUE)
  ans$order <- unlist(c)

  # Return the answer
  ans
}
