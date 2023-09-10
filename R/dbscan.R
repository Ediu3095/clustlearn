#' @title Density Based Spatial Clustering of Applications with Noise (DBSCAN)
#'
#' @description Perform DBSCAN clustering on a data matrix.
#'
#' @param data a set of observations, presented as a matrix-like object where
#' every row is a new observation.
#' @param eps how close two observations have to be to be considered neighbors.
#' @param min_pts the minimum amount of neighbors for a region to be considered
#' dense.
#' @param details a Boolean determining whether intermediate logs explaining how
#' the algorithm works should be printed or not.
#' @param waiting a Boolean determining whether the intermediate logs should be
#' printed in chunks waiting for user input before printing the next or not.
#' @param ... additional arguments passed to [proxy::dist()].
#'
#' @details The data given by \code{data} is clustered by the DBSCAN method,
#' which aims to partition the points into clusters such that the points in a
#' cluster are close to each other and the points in different clusters are far
#' away from each other. The clusters are defined as dense regions of points
#' separated by regions of low density.
#'
#' The DBSCAN method follows a 2 step process:
#'
#' \enumerate{
#'  \item For each point, the neighborhood of radius \code{eps} is computed. If
#'  the neighborhood contains at least \code{min_pts} points, then the point is
#'  considered a \strong{core point}. Otherwise, the point is considered an
#'  \strong{outlier}.
#'  \item For each core point, if the core point is not already assigned to a
#'  cluster, a new cluster is created and the core point is assigned to it.
#'  Then, the neighborhood of the core point is explored. If a point in the
#'  neighborhood is a core point, then the neighborhood of that point is also
#'  explored. This process is repeated until all points in the neighborhood have
#'  been explored. If a point in the neighborhood is not already assigned to a
#'  cluster, then it is assigned to the cluster of the core point.
#' }
#'
#' Whatever points are not assigned to a cluster are considered outliers.
#'
#' @return A [clustlearn::dbscan()] object. It is a list with the following
#' components:
#' \tabular{ll}{
#'  \code{cluster} \tab a vector of integers (from 0 to \code{max(cl$cluster)})
#'  indicating the cluster to which each point belongs. Points in cluster number
#'  0 are considered outliers. \cr
#'  \code{eps} \tab the value of \code{eps} used. \cr
#'  \code{min_pts} \tab the value of \code{min_pts} used. \cr
#'  \code{size} \tab a vector with the number of data points belonging to each
#'  cluster (where the first element is the number of outliers). \cr
#' }
#'
#' @examples
#' ### Helper function
#' test <- function(db, eps) {
#'   print(cl <- clustlearn::dbscan(db, eps))
#'   out <- cl$cluster == 0
#'   plot(db[!out, ], col = cl$cluster[!out], pch = 20, asp = 1)
#'   points(db[out, ], col = max(cl$cluster) + 1, pch = 4, lwd = 2)
#' }
#'
#' ### Example 1
#' test(clustlearn::db1, 0.3)
#'
#' ### Example 2
#' test(clustlearn::db2, 0.3)
#'
#' ### Example 3
#' test(clustlearn::db3, 0.25)
#'
#' ### Example 4
#' test(clustlearn::db4, 0.2)
#'
#' ### Example 5
#' test(clustlearn::db5, 0.3)
#'
#' ### Example 6
#' test(clustlearn::db6, 0.3)
#'
#' ### Example 7 (with explanations, no plots)
#'   cl <- clustlearn::dbscan(
#'   clustlearn::db5[1:20, ],
#'   0.3,
#'   details = TRUE,
#'   waiting = FALSE
#' )
#'
#' @author Eduardo Ruiz Sabajanes, \email{eduardo.ruizs@@edu.uah.es}
#'
#' @importFrom proxy dist
#' @export
dbscan <- function(
  data,
  eps,
  min_pts = 4,
  details = FALSE,
  waiting = TRUE,
  ...
) {
  if (details) {
    hline()
    console.log("EXPLANATION:")
    console.log("")
    console.log("The data given by data is clustered by the DBSCAN method, which aims to partition the points into clusters such that the points in a cluster are close to each other and the points in different clusters are far away from each other. The clusters are defined as dense regions of points separated by regions of low density.")
    console.log("")
    console.log("The DBSCAN method follows a 2 step process:")
    console.log("")
    console.log("    1. For each point, the neighborhood of radius eps is computed. If the neighborhood contains at least min_pts points, then the point is considered a core point. Otherwise, the point is considered an outlier.")
    console.log("    2. For each core point, if the core point is not already assigned to a cluster, a new cluster is created and the core point is assigned to it. Then, the neighborhood of the core point is explored. If a point in the neighborhood is a core point, then the neighborhood of that point is also explored. This process is repeated until all points in the neighborhood have been explored. If a point in the neighborhood is not already assigned to a cluster, then it is assigned to the cluster of the core point.")
    console.log("")
    console.log("Whatever points are not assigned to a cluster are considered outliers.")
    console.log("")

    if (waiting) {
      invisible(readline(prompt = "Press [enter] to continue"))
      console.log("")
    }
  }

  # Precompute neighbors
  distances <- as.matrix(proxy::dist(data, ...))
  neighbors <- distances <= eps

  if (details) {
    hline()
    console.log("STEP 1:")
    console.log("")
    console.log("The pairwise distances between observations are precomputed in order to later determine which of them are core observations. The distance matrix is:")
    cat("Distances:\n")
    print(round(distances, 3))
    console.log("")

    if (waiting) {
      invisible(readline(prompt = "Press [enter] to continue"))
      console.log("")
    }
  }

  # Initialize clusters
  cluster_id <- new.env()
  cluster_id$data <- data
  cluster_id$current <- 1
  cluster_id$of <- rep(-1, nrow(data))

  if (details) {
    hline()
    console.log("STEP 2:")
    console.log("")
    console.log("Every observation is labeled as UNVISITED. We are now going to loop over every observation and, if it is not already assigned to a cluster, we will try to expand a new cluster around it...")
    console.log("")

    if (waiting) {
      invisible(readline(prompt = "Press [enter] to continue"))
      console.log("")
    }
  }

  # Each loop finds a new cluster around a core point
  for (idx in seq_len(nrow(data))) {
    if (cluster_id$of[idx] != -1)
      next
    if (expand_cluster(neighbors, cluster_id, idx, min_pts, details, waiting))
      cluster_id$current <- cluster_id$current + 1
  }
  
  if (details) {
    hline()
    console.log("RESULTS:")
    console.log("")
    console.log("Having gone through every observation the following clusters have been found:")
    cat("CLUSTER #0 (NOISE):\n")
    print(cluster_id$data[cluster_id$of == 0, ])
    for (i in seq_len(max(cluster_id$of))) {
      console.log("")
      cat(paste0("CLUSTER #", i, ":\n"))
      print(cluster_id$data[cluster_id$of == i, ])
    }
    console.log("")

    if (waiting) {
      invisible(readline(prompt = "Press [enter] to continue"))
      console.log("")
    }
    hline()
  }

  # Compute cluster sizes
  size <- as.integer(table(cluster_id$of))

  # Return a dbscan object
  structure(
    list(
      cluster = cluster_id$of,
      eps = eps,
      min_pts = min_pts,
      size = size
    ),
    class = "dbscan"
  )
}

expand_cluster <- function(
  neighbors,
  cluster_id,
  point,
  min_pts,
  details,
  waiting
) {
  # Get the point's neighbors (including itself)
  seeds <- region_query(neighbors, point)

  if (length(seeds) < min_pts) {
    if (details) {
      hline()
      console.log("NOISE:")
      console.log("")
      console.log("An UNVISITED observation is labeled as NOISE:")
      cat(paste0("Observation #", point, " [UNVISITED -> NOISE]\n"))
      # print(cluster_id$data[point, ])
      console.log("")

      # console.log("With neighborhood:")
      # for (i in seeds) {
      #   tmp <- if (cluster_id$of[i] == -1) {
      #     "UNVISITED"
      #   } else if (cluster_id$of[i] == 0) {
      #     "NOISE"
      #   } else {
      #     paste0("CLUSTER #", cluster_id$of[i])
      #   }
      #   ifelse(cluster_id$of[i] == -1, "UNVISITED", "NOISE")
      #   cat(paste0("Observation #", i, " [", tmp, "]\n"))
      #   print(cluster_id$data[i, ])
      # }
      # console.log("")

      if (waiting) {
        invisible(readline(prompt = "Press [enter] to continue"))
        console.log("")
      }
    }

    # If it is not a core point, it is noise
    cluster_id$of[point] <- 0
    FALSE
  } else {
    if (details) {
      hline()
      console.log(paste0("CLUSTER #", cluster_id$current, ":"))
      console.log("")
      console.log("A new cluster is going to be expanded around an UNVISITED core observation:")
      cat(paste0("Observation #", point, " [UNVISITED -> CLUSTER #", cluster_id$current, "]\n"))
      # print(cluster_id$data[point, ])
      console.log("")

      # console.log("With neighborhood:")
      # for (i in seeds) {
      #   tmp <- if (cluster_id$of[i] == -1) {
      #     "UNVISITED"
      #   } else if (cluster_id$of[i] == 0) {
      #     "NOISE"
      #   } else {
      #     paste0("CLUSTER #", cluster_id$of[i])
      #   }
      #   ifelse(cluster_id$of[i] == -1, "UNVISITED", "NOISE")
      #   cat(paste0("Observation #", i, " [", tmp, "]\n"))
      #   print(cluster_id$data[i, ])
      # }
      # console.log("")

      if (waiting) {
        invisible(readline(prompt = "Press [enter] to continue"))
        console.log("")
      }
    }

    # Otherwise, we can expand the cluster
    frontier <- setdiff(seeds, point)

    if (details) {
      console.log("The cluster is also expanded around the neighbors of the core observation:")
      for (i in frontier) {
        tmp1 <- ifelse(cluster_id$of[i] == -1, "UNVISITED", "NOISE")
        tmp2 <- paste0("CLUSTER #", cluster_id$current)
        cat(paste0("Observation #", i, " [", tmp1, " -> ", tmp2, "]\n"))
        # print(cluster_id$data[i, ])
      }
      console.log("")

      console.log("All of these observations are added to the cluster.")
      console.log("")

      if (waiting) {
        invisible(readline(prompt = "Press [enter] to continue"))
        console.log("")
      }
    }

    cluster_id$of[seeds] <- cluster_id$current

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

        if (details) {
          console.log("***")
          console.log("")
          console.log("The following core observation is expanded:")
          cat(paste0("Observation #", current_point, " [CLUSTER #", cluster_id$current, "]\n"))
          # print(cluster_id$data[current_point, ])
          console.log("")

          console.log("It's neighborhood is:")
          for (i in result) {
            tmp <- if (cluster_id$of[i] == -1) {
              "UNVISITED"
            } else if (cluster_id$of[i] == 0) {
              "NOISE"
            } else {
              paste0("CLUSTER #", cluster_id$of[i])
            }
            ifelse(cluster_id$of[i] == -1, "UNVISITED", "NOISE")
            cat(paste0("Observation #", i, " [", tmp, "]\n"))
            # print(cluster_id$data[i, ])
          }
          console.log("")

          if (length(result[not_visited | noise]) > 0) {
            console.log("Upon doing it, the following observations are added to the cluster:")
            for (i in result[not_visited | noise]) {
              tmp1 <- ifelse(cluster_id$of[i] == -1, "UNVISITED", "NOISE")
              tmp2 <- paste0("CLUSTER #", cluster_id$current)
              cat(paste0("Observation #", i, " [", tmp1, " -> ", tmp2, "]\n"))
              # print(cluster_id$data[i, ])
            }
          } else {
            console.log("Upon doing it, no observations are added to the cluster...")
          }
          console.log("")

          if (length(result[not_visited]) > 0) {
            console.log("Additionally, these observations are also expanded:")
            for (i in result[not_visited]) {
              cat(paste0("Observation #", i, " [CLUSTER #", cluster_id$current, "]\n"))
              # print(cluster_id$data[i, ])
            }
          } else {
            console.log("Additionally, no other observations are expanded...")
          }
          console.log("")

          if (waiting) {
            invisible(readline(prompt = "Press [enter] to continue"))
            console.log("")
          }
        }

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
