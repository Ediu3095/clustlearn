#' @title Agglomerative Hierarchical Clustering
#'
#' @description Perform a hierarchical agglomerative cluster analysis on a set
#' of observations
#'
#' @param data a set of observations, presented as a matrix-like object where
#' every row is a new observation.
#' @param proximity the proximity definition to be used. This should be one
#' of \code{"single"} (minimum/single linkage), \code{"complete"} (maximum/
#' complete linkage), \code{"average"} (average linkage).
#' @param details a Boolean determining whether intermediate logs explaining how
#' the algorithm works should be printed or not.
#' @param waiting a Boolean determining whether the intermediate logs should be
#' printed in chunks waiting for user input before printing the next or not.
#' @param ... additional arguments passed to [proxy::dist()].
#'
#' @details This function performs a hierarchical cluster analysis for the
#' \eqn{n} objects being clustered. The definition of a set of clusters using
#' this method follows a \eqn{n} step process, which repeats until a single
#' cluster remains:
#'
#' \enumerate{
#'  \item Initially, each object is assigned to its own cluster. The matrix
#'  of distances between clusters is computed.
#'  \item The two clusters with closest proximity will be joined together and
#'  the proximity matrix updated. This is done according to the specified
#'  \code{proximity}. This step is repeated until a single cluster remains.
#' }
#'
#' The definitions of \code{proximity} considered by this function are:
#'
#' \describe{
#'  \item{\code{single}}{\eqn{\min\left\{d(x,y):x\in A,y\in B\right\}}. Defines
#'  the proximity between two clusters as the distance between the closest
#'  objects among the two clusters. It produces clusters where each object is
#'  closest to at least one other object in the same cluster. It is known as
#'  \strong{SLINK}, \strong{single-link} and \strong{minimum-link}.}
#'  \item{\code{complete}}{\eqn{\max\left\{d(x,y):x\in A,y\in B\right\}}.
#'  Defines the proximity between two clusters as the distance between the
#'  furthest objects among the two clusters. It is known as \strong{CLINK},
#'  \strong{complete-link} and \strong{maximum-link}.}
#'  \item{\code{average}}{\eqn{\frac{1}{\left|A\right|\cdot\left|B\right|}
#'  \sum_{x\in A}\sum_{y\in B} d(x,y)}. Defines the proximity between two
#'  clusters as the average distance between every pair of objects, one from
#'  each cluster. It is also known as \strong{UPGMA} or \strong{average-link}.}
#' }
#'
#' @return An [stats::hclust()] object which describes the tree produced by the
#' clustering process.
#'
#' @examples
#' ### !! This algorithm is very slow, so we'll only test it on some datasets !!
#'
#' ### Helper function
#' test <- function(db, k, prox) {
#'   print(cl <- clustlearn::agglomerative_clustering(db, prox))
#'   oldpar <- par(mfrow = c(1, 2))
#'   plot(db, col = cutree(cl, k), asp = 1, pch = 20)
#'   h <- rev(cl$height)[50]
#'   clu <- as.hclust(cut(as.dendrogram(cl), h = h)$upper)
#'   ctr <- unique(cutree(cl, k)[cl$order])
#'   plot(clu, labels = FALSE, hang = -1, xlab = "Cluster", sub = "", main = "")
#'   rect.hclust(clu, k = k, border = ctr)
#'   par(oldpar)
#' }
#'
#' ### Example 1
#' test(clustlearn::db1, 2, "single")
#'
#' ### Example 2
#' # test(clustlearn::db2, 2, "sing") # same as "single"
#'
#' ### Example 3
#' test(clustlearn::db3, 4, "a") # same as "average"
#'
#' ### Example 4
#' test(clustlearn::db4, 6, "s") # same as "single"
#'
#' ### Example 5
#' test(clustlearn::db5, 3, "complete")
#'
#' ### Example 6
#' # test(clustlearn::db6, 3, "c") # same as "complete"
#'
#' ### Example 7 (with explanations, no plots)
#'   cl <- clustlearn::agglomerative_clustering(
#'   clustlearn::db5[1:6, ],
#'   'single',
#'   details = TRUE,
#'   waiting = FALSE
#' )
#'
#' @author Eduardo Ruiz Sabajanes, \email{eduardo.ruizs@@edu.uah.es}
#'
#' @importFrom proxy dist
#' @importFrom stats as.dist
#' @export
agglomerative_clustering <- function(
  data,
  proximity = "single",
  details = FALSE,
  waiting = TRUE,
  ...
) {
  # Function needed to calculate the avg distance between two clusters
  avg <- function(m1, m2) function(d1, d2) (d1 * m1 + d2 * m2) / (m1 + m2)

  # Figure out the proximity definition
  proximity <- grep(
    tolower(proximity),
    c("single", "complete", "average"),
    value = TRUE,
    fixed = TRUE
  )

  # Exactly one proximity definition should be found
  if (length(proximity) != 1)
    stop("Invalid proximity method")

  if (details) {
    hline()
    console.log("EXPLANATION:")
    console.log("")
    console.log("The Agglomerative Hierarchical Clustering algorithm defines a clustering hierarchy for a dataset following a `n` step process, which repeats until a single cluster remains:")
    console.log("")
    console.log("    1. Initially, each object is assigned to its own cluster. The matrix of distances between clusters is computed.")
    console.log("    2. The two clusters with closest proximity will be joined together and the proximity matrix updated. This is done according to the specified proximity. This step is repeated until a single cluster remains.")
    console.log("")
    console.log("The definitions of proximity considered by this function are:")
    console.log("")
    console.log("    1. `single`. Defines the proximity between two clusters as the distance between the closest objects among the two clusters. It produces clusters where each object is closest to at least one other object in the same cluster. It is known as SLINK, single-link or minimum-link.")
    console.log("    2. `complete`. Defines the proximity between two clusters as the distance between the furthest objects among the two clusters. It is known as CLINK, complete-link or maximum-link.")
    console.log("    3. `average`. Defines the proximity between two clusters as the average distance between every pair of objects, one from each cluster. It is also known as UPGMA or average-link.")
    console.log("")

    if (waiting) {
      invisible(readline(prompt = "Press [enter] to continue"))
      console.log("")
    }
  }

  # Prepare the data structure which will hold the final answer
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

  # Create a list with the initial clusters
  cl <- lapply(
    seq_len(nrow(data)),
    function(data) {
      structure(
        data,
        label = -data,
        members = 1
      )
    }
  )
  tmp <- sapply(cl, function(x) attr(x, "label"))

  if (details) {
    hline()
    console.log("STEP 1:")
    console.log("")
    console.log("Initially, each object is assigned to its own cluster. This leaves us with the following clusters:")
    for (i in seq_len(length(cl))) {
      cat(paste0("CLUSTER #", attr(cl[[i]], "label"), " (size: ", attr(cl[[i]], "members"), ")", "\n"))
      print(data[cl[[i]], , drop = FALSE])
    }
    console.log("")

    if (waiting) {
      invisible(readline(prompt = "Press [enter] to continue"))
      console.log("")
    }
  }

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
  dimnames(d) <- list(tmp, tmp)

  if (details) {
    console.log("The matrix of distances between clusters is computed:")
    cat("Distances:\n")
    print(as.dist(round(d, 3)))
    console.log("")

    if (waiting) {
      invisible(readline(prompt = "Press [enter] to continue"))
      console.log("")
    }
  }

  for (i in seq_len(length(cl) - 1)) {
    # Look for the minimum distance between two clusters
    md <- which.min(d) - 1
    md <- sort(c(md %% nrow(d), md %/% nrow(d)) + 1)

    # Join the clusters into a new one
    c1 <- cl[[md[1]]]
    m1 <- attr(c1, "members")
    c2 <- cl[[md[2]]]
    m2 <- attr(c2, "members")
    c3 <- structure(
      list(c1, c2),
      label = i,
      members = m1 + m2
    )

    if (details) {
      hline()
      console.log(paste0("STEP ", i + 1, ":"))
      console.log("")
      console.log("The two clusters with closest proximity are identified:")
      cat("Clusters:\n")
      cat(paste0("CLUSTER #", attr(c1, "label"), " (size: ", m1, ")", "\n"))
      cat(paste0("CLUSTER #", attr(c2, "label"), " (size: ", m2, ")", "\n"))
      cat("Proximity:\n")
      print(d[md[1], md[2]])
      console.log("")

      if (waiting) {
        invisible(readline(prompt = "Press [enter] to continue"))
        console.log("")
      }
    }

    if (details) {
      console.log("They are merged into a new cluster:")
      cat(paste0("CLUSTER #", attr(c3, "label"), " (size: ", attr(c3, "members"), ") [", "CLUSTER #", attr(c1, "label"), " + CLUSTER #", attr(c2, "label"), "]\n"))
      console.log("")

      if (waiting) {
        invisible(readline(prompt = "Press [enter] to continue"))
        console.log("")
      }
    }

    # Add the merged clusters to the answer
    ans$merge <- c(ans$merge, attr(c1, "label"), attr(c2, "label"))
    ans$height <- c(ans$height, d[md[1], md[2]])
    cl <- c(cl, list(c3))
    cl <- cl[-md]
    tmp <- sapply(cl, function(x) attr(x, "label"))

    # Recompute the distances (proximity)
    d1 <- d[, md[1]]
    d2 <- d[, md[2]]
    d3 <- mapply(
      switch(
        proximity,
        single = min,
        complete = max,
        average = avg(m1, m2)
      ),
      d1,
      d2
    )
    d3[md] <- Inf
    d <- cbind(d, d3)
    d <- rbind(d, c(d3, Inf))
    d <- d[-md, -md, drop = FALSE]
    dimnames(d) <- list(tmp, tmp)

    if (details) {
      console.log("The proximity matrix is updated. To do so the rows/columns of the merged clusters are removed, and the rows/columns of the new cluster are added:")
      cat("Distances:\n")
      print(as.dist(round(d, 3)))
      console.log("")

      if (waiting) {
        invisible(readline(prompt = "Press [enter] to continue"))
        console.log("")
      }
    }
  }

  # Compute the merge and order of the hclust
  ans$merge <- matrix(ans$merge, ncol = 2, byrow = TRUE)
  ans$order <- unlist(cl)

  if (details) {
    hline()
    console.log("RESULTS:")
    console.log("")
    console.log("Since all clusters have been merged together, the final clustering hierarchy is:")
    console.log("(Check the plot for the dendrogram representation of the hierarchy)")
    plot(ans, hang = -1)
    console.log("")

    if (waiting) {
      invisible(readline(prompt = "Press [enter] to continue"))
      console.log("")
    }

    hline()
  }

  # Return the answer
  ans
}
