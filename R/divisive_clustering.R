#' @title Divisive Hierarchical Clustering
#'
#' @description Perform a hierarchical Divisive cluster analysis on a set of
#' observations
#'
#' @param data a set of observations, presented as a matrix-like object where
#' every row is a new observation.
#' @param ... additional arguments passed to [clustlearn::kmeans()].
#'
#' @details This function performs a hierarchical cluster analysis for the
#' \eqn{n} objects being clustered. The definition of a set of clusters using
#' this method follows a \eqn{n} step process, which repeats until \eqn{n}
#' clusters remain:
#'
#' \enumerate{
#'  \item Initially, each object is assigned to the same cluster. The sum of
#'  squares of the distances between objects and their cluster center is
#'  computed.
#'  \item The cluster with the highest sum of squares is split into two using
#'  the k-means algorithm. This step is repeated until \eqn{n} clusters remain.
#' }
#'
#' @return An [stats::hclust()] object which describes the tree produced by the
#' clustering process.
#'
#' @examples
#' ### Helper function
#' test <- function(db, k) {
#'   print(cl <- clustlearn::divisive_clustering(db))
#'   par(mfrow = c(1, 2))
#'   plot(db, col = cutree(cl, k), asp = 1, pch = 20)
#'   h <- rev(cl$height)[50]
#'   clu <- as.hclust(cut(as.dendrogram(cl), h = h)$upper)
#'   ctr <- unique(cutree(cl, k)[cl$order])
#'   plot(clu, labels = FALSE, hang = -1, xlab = "Cluster", sub = "", main = "")
#'   rect.hclust(clu, k = k, border = ctr)
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
#' @export
divisive_clustering <- function(data, ...) {
  # Prepare the data structure which will hold the answer
  ans <- structure(
    list(
      merge = numeric(0),
      height = NULL,
      order = NULL,
      labels = rownames(data),
      method = "kmeans",
      call = NULL,
      dist.method = "Euclidean"
    ),
    class = "hclust"
  )

  # Wrap the data with additional information we'll need
  data_center <- apply(data, 2, mean)
  totss <- sum(apply(data, 1, function(x) x - data_center)^2)
  wrapped_data <- list(
    data = data,
    label = nrow(data) - 1,
    ss = totss,
    elems = -seq_len(nrow(data))
  )

  # Build a list with all clusters
  clusters <- list(wrapped_data)

  # Until there are no clusters with sum of squares greater than 0
  while (any(sapply(clusters, function(x) length(x$elems)) > 1)) {
    # We'll operate on the first cluster and order them all afterwards
    target <- clusters[[1]]

    # Split the target cluster into two using the k-means approach
    km <- clustlearn::kmeans(target$data, 2, ...)
    if (any(km$size == 0))
      km$cluster[1] <- setdiff(1:2, km$cluster[1])

    # Create clusters for each split
    lhs <- list(
      data = target$data[km$cluster == 1, , drop = FALSE],
      label = NULL,
      ss = km$withinss[1],
      elems = target$elems[km$cluster == 1]
    )
    lhs$label <- if (length(lhs$elems) == 1) {
      lhs$elems
    } else {
      target$label - 1
    }
    rhs <- list(
      data = target$data[km$cluster == 2, , drop = FALSE],
      label = NULL,
      ss = km$withinss[2],
      elems = target$elems[km$cluster == 2]
    )
    rhs$label <- if (length(rhs$elems) == 1) {
      rhs$elems
    } else {
      target$label - length(lhs$elems)
    }

    # Update the answer
    ans$merge <- c(lhs$label, rhs$label, ans$merge)
    ans$height <- c(km$totss, ans$height)

    # Replace the target cluster with the two new ones
    clusters <- c(
      if (length(lhs$elems) > 1 && length(rhs$elems) > 1) {
        list(lhs, rhs)
      } else if (length(lhs$elems) > 1) {
        list(lhs)
      } else if (length(rhs$elems) > 1) {
        list(rhs)
      } else {
        list()
      },
      clusters[-1]
    )
  }

  # Compute the merge and order of the hclust
  ans$merge <- matrix(ans$merge, ncol = 2, byrow = TRUE)
  ans$height <- sqrt(ans$height)

  # Order the rows in merge by height
  tmp <- order_merge_by_height(ans$merge, ans$height)
  ans$merge <- tmp$merge
  ans$height <- tmp$height

  # Figure the order of the elements from the merge to build a proper dendrogram
  ans$order <- merge2order(ans$merge)

  # Return the answer
  ans
}

order_merge_by_height <- function(merge, height) {
  ord1 <- order(height)
  ord2 <- order(ord1)
  for (i in seq_len(nrow(merge))) {
    for (j in seq_len(ncol(merge))) {
      val <- merge[i, j]
      merge[i, j] <- if (val < 0) val else ord2[val]
    }
  }
  merge <- merge[ord1, , drop = FALSE]
  height <- height[ord1]

  list(merge = merge, height = height)
}

merge2order <- function(merge) {
  order <- if (nrow(merge) > 0) merge[nrow(merge), ] else -1
  while (any(order > 0)) {
    target <- which.max(order)
    order <- c(
      order[seq_along(order) < target],
      merge[order[target], ],
      order[seq_along(order) > target]
    )
  }
  abs(order)
}
