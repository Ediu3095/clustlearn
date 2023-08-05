#' @title Divisive Hierarchical Clustering
#'
#' @description Hierarchical divisive cluster analysis on a set of observations
#'
#' @param data a set of observations, presented as a matrix-like object where
#' every row is a new observation.
#' @param ... additional arguments passed to [kmeans()].
#'
#' @details
#'
#' @author Eduardo Ruiz Sabajanes, \email{eduardo.ruizs@@edu.uah.es}
#'
#' @return An [hclust()] object which describes the tree produced by the
#' clustering process.
#'
#' @examples
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
  label_max <- wrapped_data$label
  while (any(sapply(clusters, function(x) length(x$elems)) > 1)) {
    # We'll operate on the cluster with greatest sum of squares
    target <- which.max(sapply(clusters, function(x) x$ss))

    # Split the target cluster into two using the k-means approach
    kmeans_split <- clustlearn::kmeans(clusters[[target]]$data, 2, max_iterations)
    if (length(unique(kmeans_split$cluster)) < 2)
      kmeans_split$cluster[1] <- 2

    # Create clusters for each split
    lhs <- list(
      data = clusters[[target]]$data[kmeans_split$cluster == 1, , drop = FALSE],
      label = NULL,
      ss = kmeans_split$withinss[1],
      elems = clusters[[target]]$elems[kmeans_split$cluster == 1]
    )
    lhs$label <- if (nrow(lhs$data) == 1) {
      lhs$elems
    } else {
      label_max <- label_max - 1
    }
    rhs <- list(
      data = clusters[[target]]$data[kmeans_split$cluster == 2, , drop = FALSE],
      label = NULL,
      ss = kmeans_split$withinss[2],
      elems = clusters[[target]]$elems[kmeans_split$cluster == 2]
    )
    rhs$label <- if (nrow(rhs$data) == 1) {
      rhs$elems
    } else {
      label_max <- label_max - 1
    }

    # Update the answer
    ans$merge <- c(lhs$label, rhs$label, ans$merge)
    ans$height <- c(kmeans_split$totss, ans$height)

    # Replace the target cluster with the two new ones
    clusters <- c(
      clusters[seq_along(clusters) < target],
      if (length(lhs$elems) > 1 && length(rhs$elems) > 1) {
        list(lhs, rhs)
      } else if (length(lhs$elems) > 1) {
        list(lhs)
      } else if (length(rhs$elems) > 1) {
        list(rhs)
      } else {
        list()
      },
      clusters[seq_along(clusters) > target]
    )
  }

  # Compute the merge and order of the hclust
  ans$merge <- matrix(ans$merge, ncol = 2, byrow = TRUE)
  # ans$height <- seq_len(nrow(data)-1)
  ans$height <- sqrt(ans$height)
  # ans$order <- abs(sapply(clusters, function(x) x$label))
  ans$order <- merge2order(ans$merge)

  # Return the answer
  ans
}

merge2order <- function(merge) {
  order <- if (nrow(merge) > 0) { merge[nrow(merge), ] } else { -1 }
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
