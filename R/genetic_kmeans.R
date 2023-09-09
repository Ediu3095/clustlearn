#' @title K-Means Clustering
#'
#' @description Perform K-Means clustering on a data matrix.
#'
#' @param data a set of observations, presented as a matrix-like object where
#' every row is a new observation.
#' @param k the number of clusters.
#' @param population_size the number of individuals in the population.
#' @param mut_probability the probability of a mutation occurring.
#' @param max_generations the maximum number of iterations allowed.
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
#'    This is determined through the use of the \code{k} parameter.
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
#'  \item{\code{random}}{A set of \code{k} observations is chosen at
#'  random from the data as the initial centers.}
#'  \item{\code{kmeans++}}{The \code{k} observations are chosen using the
#'  \strong{kmeans++} algorithm. This algorithm chooses the first center at
#'  random and then chooses the next center from the remaining observations with
#'  probability proportional to the square distance to the closest center. This
#'  process is repeated until \code{k} centers are chosen.}
#' }
#'
#' @return A [stats::kmeans()] object.
#'
#' @examples
#' ### Voronoi tesselation
#' voronoi <- require(deldir)
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
#' @author Eduardo Ruiz Sabajanes, \email{eduardo.ruizs@@edu.uah.es}
#'
#' @importFrom proxy dist
#' @importFrom stats runif sd
#' @export
genetic_kmeans <- function(
  data,
  k,
  population_size = 10,
  mut_probability = .5,
  max_generations = 10,
  ...
) {
  # Initialize the population
  population <- gka_initialization(nrow(data), population_size, k)

  # Compute the fitness of the initial population
  centers <- gka_centers(data, k, population)
  twcv <- gka_twcv(data, k, population, centers)
  fitness <- gka_fitness(twcv)

  # Choose the best individual as the initial solution
  idx <- which.max(fitness)
  solution <- structure(
    list(
      cluster = population[idx, ],
      centers = centers[idx, , ],
      totss = NULL,
      withinss = NULL,
      tot.withinss = twcv[idx],
      betweenss = NULL,
      size = NULL,
      iter = NULL,
      ifault = 0
    ),
    class = "kmeans"
  )

  # Run the algorithm for a maximum of max_generations generations or until the
  # solution converges to a local optimum (i.e. the fitness of the best
  # individual does not improve)
  iter <- 0
  while (iter < max_generations) {
    iter <- iter + 1

    # Selection
    idx <- gka_selection(population_size, fitness)
    chromosome <- population[idx, ]
    centers <- centers[idx, , ]

    # Mutation
    prob <- gka_allele_mutation(data, k, centers, ...)
    for (i in seq_len(population_size)) {
      population[i, ] <- gka_mutation(chromosome, prob, k, mut_probability)
    }
    population <- gka_chromosome_fix(population, k)

    # Crossover
    centers <- gka_centers(data, k, population)
    for (i in seq_len(population_size)) {
      population[i, ] <- gka_crossover(data, centers[i, , ])
    }
    population <- gka_chromosome_fix(population, k)

    # Compute the fitness of the new population
    centers <- gka_centers(data, k, population)
    twcv <- gka_twcv(data, k, population, centers)
    fitness <- gka_fitness(twcv)
    idx <- which.max(fitness)

    # Update the solution
    if (twcv[idx] < solution$tot.withinss) {
      idx <- which.max(fitness)
      solution$cluster <- population[idx, ]
      solution$centers <- centers[idx, , ]
      solution$tot.withinss <- twcv[idx]
    }
  }

  # Name the centroids
  row.names(solution$centers) <- seq_len(k)

  # Total sum of squares
  center <- apply(data, 2, mean)
  totss <- sum(apply(data, 1, function(x) x - center) ^ 2)

  # Total within-cluster sum of squares
  withinss <- sapply(
    seq_len(k),
    function(cluster) {
      ccenter <- solution$centers[cluster, ]
      cdata <- data[solution$cluster == cluster, , drop = FALSE]
      sum(apply(cdata, 1, function(x) x - ccenter)^2)
    }
  )

  # Total within-cluster sum of squares
  tot_withinss <- sum(withinss)

  # The between-cluster sum of squares
  betweenss <- totss - solution$tot.withinss

  # Find the size of each cluster
  tmp <- factor(solution$cluster, levels = seq_len(k))
  size <- as.integer(table(tmp))

  # Update the solution
  solution$totss <- totss
  solution$withinss <- withinss
  solution$betweenss <- betweenss
  solution$size <- size
  solution$iter <- iter

  solution
}

# @title Initialization method
#
# @param n the number of observations in the data.
# @param p the number of individuals in the population.
# @param k the number of clusters.
#
# @return a matrix of size \code{p} by \code{n} with the cluster assignments
# for each observation.
#
# @author Eduardo Ruiz Sabajanes, \email{eduardo.ruizs@@edu.uah.es}
gka_initialization <- function(n, p, k) {
  population <- matrix(0, nrow = p, ncol = n)

  # For each member of the population
  for (i in seq_len(p)) {
    # Assign p observations to each cluster, where p is the floor of n / k
    quotient <- n %/% k
    pigeons <- sample(seq_len(n), quotient * k)
    for (j in seq_len(k)) {
      from <- (j - 1) * quotient + 1
      to   <- from    + quotient - 1
      population[i, pigeons[seq(from, to)]] <- j
    }

    # Assign the remaining observations to a random cluster
    remainder <- n %% k
    population[i, -pigeons] <- sample(seq_len(k), remainder)
  }

  population
}

# @title Centroid computation
#
# @param data a set of observations, presented as a matrix-like object where
# every row is a new observation. The matrix is of size \code{n} by \code{m}.
# @param k the number of clusters.
# @param population a matrix of size \code{p} by \code{n} with the cluster
# assignments for each observation.
#
# @return a 3D array of size \code{p} by \code{k} by \code{m}.
#
# @author Eduardo Ruiz Sabajanes, \email{eduardo.ruizs@@edu.uah.es}
gka_centers <- function(data, k, population) {
  centers <- array(0, dim = c(nrow(population), k, ncol(data)))
  for (i in seq_len(nrow(population))) {
    for (j in seq_len(k)) {
      centers[i, j, ] <- colMeans(data[population[i, ] == j, ])
    }
  }
  centers
}

# @title Total Within Cluster Variation (TWCV) computation
#
# @param data a set of observations, presented as a matrix-like object where
# every row is a new observation. The matrix is of size \code{n} by \code{m}.
# @param k the number of clusters.
# @param population a matrix of size \code{p} by \code{n} with the cluster
# assignments for each observation.
# @param centers a 3D array of size \code{p} by \code{k} by \code{m} with the
# cluster centers for each individual in the population.
#
# @return a vector of size \code{p} with the total within cluster variation of
# each individual in the population.
#
# @author Eduardo Ruiz Sabajanes, \email{eduardo.ruizs@@edu.uah.es}
gka_twcv <- function(data, k, population, centers) {
  twcv <- numeric(nrow(population))
  for (i in seq_len(nrow(population))) {
    twcv[i] <- sum(
      sapply(
        seq_len(k),
        function(j) sum((t(data[population[i, ] == j, ]) - centers[i, j, ]) ^ 2)
      )
    )
  }
  twcv
}

# @title Fitness function
#
# @param twcv a vector of size \code{p} with the total within cluster
# variation of each individual in the population.
#
# @return a vector of size \code{p} with the fitness of each individual in the
# population.
#
# @author Eduardo Ruiz Sabajanes, \email{eduardo.ruizs@@edu.uah.es}
gka_fitness <- function(twcv) {
  f <- -twcv
  g <- f - (mean(f) - 2 * sd(f))
  g[g <= 0] <- 1e-6
  g
}

# @title Selection method
#
# @param p the number of individuals in the population.
# @param fitness a vector of size \code{p} with the fitness of each individual
# in the population.
#
# @return the index of the individual selected for reproduction.
#
# @author Eduardo Ruiz Sabajanes, \email{eduardo.ruizs@@edu.uah.es}
gka_selection <- function(p, fitness) {
  sample(seq_len(p), 1, prob = fitness / sum(fitness))
}

# @title Allele mutation probability computation
#
# @param data a set of observations, presented as a matrix-like object where
# every row is a new observation. The matrix is of size \code{n} by \code{m}.
# @param k the number of clusters.
# @param centers a matrix of size \code{k} by \code{m} with the cluster centers
# @param ... additional arguments passed to [proxy::dist()].
#
# @return a matrix of size \code{n} by \code{k} with the probability of each
# allele mutating to a specific cluster.
#
# @author Eduardo Ruiz Sabajanes, \email{eduardo.ruizs@@edu.uah.es}
#
# @importFrom proxy dist
gka_allele_mutation <- function(data, k, centers, ...) {
  d <- proxy::dist(data, centers, ...)
  dmax <- apply(d, 1, max)
  tmp1 <- matrix(0, nrow = nrow(data), ncol = k)
  for (j in seq_len(k)) {
    tmp1[, j] <- 1.5 * dmax - d[, j]
  }
  tmp1[tmp1 <= 0] <- 1e-6
  tmp2 <- rowSums(tmp1)
  prob <- tmp1 / tmp2
  prob
}

# @title Mutation method
#
# @param chromosome a vector of size \code{n} with the cluster assignments for
# each observation.
# @param prob a matrix of size \code{n} by \code{k} with the probability of
# each allele mutating to a specific cluster.
# @param k the number of clusters.
# @param mut_probability the probability of a mutation occurring.
#
# @return a vector of size \code{n} with the cluster assignments for each
# observation i.e. a new chromosome.
#
# @author Eduardo Ruiz Sabajanes, \email{eduardo.ruizs@@edu.uah.es}
#
# @importFrom stats runif
gka_mutation <- function(chromosome, prob, k, mut_probability) {
  new_chromosome <- chromosome
  for (i in seq_len(length(new_chromosome))) {
    if (runif(1) < mut_probability) {
      new_chromosome[i] <- sample(seq_len(k), 1, prob = prob[i, ])
    }
  }
  new_chromosome
}

# @title Crossover method i.e. K-Means Operator
#
# @description K-Means Operator (KMO) which replaces the crossover operator in
# the Genetic K-Means algorithm (GKA).
#
# @param data a set of observations, presented as a matrix-like object where
# every row is a new observation. The matrix is of size \code{n} by \code{m}.
# @param centers a matrix of size \code{k} by \code{m} with the cluster centers
# for a specific individual in the population.
#
# @return a vector of size \code{n} with the cluster assignments for each
# observation i.e. a new chromosome.
#
# @author Eduardo Ruiz Sabajanes, \email{eduardo.ruizs@@edu.uah.es}
#
# @importFrom proxy dist
gka_crossover <- function(data, centers) {
  d <- proxy::dist(data, centers)
  apply(d, 1, which.min)
}

# @title Chromosome fixing method
#
# @description This method fixes chromosomes which do not have at least one
# observation assigned to each cluster.
#
# @param population a matrix of size \code{p} by \code{n} with the cluster
# assignments for each observation.
# @param k the number of clusters.
#
# @return a matrix of size \code{p} by \code{n} with the cluster assignments
# for each observation.
#
# @author Eduardo Ruiz Sabajanes, \email{eduardo.ruizs@@edu.uah.es}
gka_chromosome_fix <- function(population, k) {
  for (i in seq_len(nrow(population))) {
    sdiff <- setdiff(seq_len(k), unique(population[i, ]))
    if (length(sdiff) <= 0)
      next

    for (cluster in sdiff) {
      n <- ncol(population)
      population[i, sample(seq_len(n), n %/% k)] <- cluster
    }
  }
  population
}
