#' @title Gaussian mixture model
#'
#' @description Perform Gaussian mixture model clustering on a data matrix.
#'
#' @param data a set of observations, presented as a matrix-like object where
#' every row is a new observation.
#' @param k the number of clusters to find.
#' @param max_iter the maximum number of iterations to perform.
#' @param ... additional arguments passed to [clustlearn::kmeans()].
#'
#' @details The data given by \code{data} is clustered by the model-based
#' algorithm that assumes every cluster follows a normal distribution, thus
#' the name "Gaussian Mixture".
#'
#' The normal distributions are parameterized by their mean vector, covariance
#' matrix and mixing proportion. Initially, the mean vector is set to the
#' cluster centers obtained by performing a k-means clustering on the data,
#' the covariance matrix is set to the covariance matrix of the data points
#' belonging to each cluster and the mixing proportion is set to the proportion
#' of data points belonging to each cluster. The algorithm then optimizes the
#' gaussian models by means of the Expectation Maximization (EM) algorithm.
#'
#' The EM algorithm is an iterative algorithm that alternates between two steps:
#'
#' \describe{
#'  \item{Expectation}{Compute the expected value of the log likelihood of the
#'  data given the current parameters.}
#'  \item{Maximization}{Update the parameters to maximize the expected value
#'  of the log likelihood.}
#' }
#'
#' The algorithm stops when the difference between the expected value of the
#' log likelihood of the data given the current parameters and the previous
#' iteration is less than 1e-6 or when the maximum number of iterations is
#' reached.
#'
#' @return A [clustlearn::gaussian_mixture()] object. It is a list with the
#' following components:
#' \tabular{ll}{
#'  \code{cluster} \tab a vector of integers (from \code{1:k}) indicating the
#'  cluster to which each point belongs. \cr
#'  \code{mu} \tab the final mean parameters. \cr
#'  \code{sigma} \tab the final covariance matrices. \cr
#'  \code{lambda} \tab the final mixing proportions. \cr
#'  \code{loglik} \tab the final log likelihood. \cr
#'  \code{all.loglik} \tab a vector of each iteration's log likelihood. \cr
#'  \code{iter} \tab the number of iterations performed. \cr
#'  \code{size} \tab a vector with the number of data points belonging to each
#'  cluster. \cr
#' }
#'
#' @examples
#' ### Helper function
#' test <- function(db, k) {
#'   print(cl <- clustlearn::gaussian_mixture(db, k, 100))
#'
#'   x <- seq(min(db[, 1]), max(db[, 1]), length.out = 100)
#'   y <- seq(min(db[, 2]), max(db[, 2]), length.out = 100)
#'
#'   plot(db, col = cl$cluster, asp = 1, pch = 20)
#'   for (i in seq_len(k)) {
#'     m <- cl$mu[i, ]
#'     s <- cl$sigma[i, , ]
#'     f <- function(x, y) cl$lambda[i] * clustlearn:::dmnorm(cbind(x, y), m, s)
#'     z <- outer(x, y, f)
#'     contour(x, y, z, col = i, add = TRUE)
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
#' @importFrom stats cov cov.wt
#' @export
gaussian_mixture <- function(data, k, max_iter = 10, ...) {
  data <- as.matrix(data)

  # Perform k-means to get initial values for mu, sigma and pi
  members <- kmeans(data, k, ...)
  mu <- members$centers
  sigma <- array(0, dim = c(k, ncol(data), ncol(data)))
  for (i in seq_len(k)) {
    sigma[i, , ] <- as.matrix(cov(data[members$cluster == i, , drop = FALSE]))
  }
  lambda <- members$size / nrow(data)

  # EM algorithm
  # Starting values of expected value of the log likelihood
  q <- c(
    sum.finite(
      sapply(
        seq_len(k),
        function(i) {
          log(lambda[i]) + log(dmnorm(data, mu[i, ], as.matrix(sigma[i, , ])))
        }
      )
    ),
    0
  )
  it <- 0
  while (abs(diff(q[1:2])) >= 1e-6 && it < max_iter) {
    # E step
    comp <- sapply(
      seq_len(k),
      function(i) lambda[i] * dmnorm(data, mu[i, ], as.matrix(sigma[i, , ]))
    )
    comp_sum <- rowSums.finite(comp)
    p <- comp / comp_sum

    # M step
    lambda <- sapply(
      seq_len(k),
      function(i) sum.finite(p[, i]) / nrow(data)
    )
    for (i in seq_len(k)) {
      mu[i, ] <- colSums.finite(p[, i] * data) / sum.finite(p[, i])
    }
    for (i in seq_len(k)) {
      tmp <- cov.wt.finite(data, wt = p[, i], center = mu[i, ])$cov
      sigma[i, , ] <- as.matrix(tmp)
    }

    # Compute new expected value of the log likelihood
    q <- c(sum(log(comp_sum)), q)
    it <- it + 1
  }

  comp <- sapply(
    seq_len(k),
    function(i) lambda[i] * dmnorm(data, mu[i, ], as.matrix(sigma[i, , ]))
  )
  cluster <- apply(comp, 1, which.max)
  size <- as.integer(table(factor(cluster, levels = seq_len(k))))

  structure(
    list(
      cluster = cluster,
      mu = mu,
      sigma = sigma,
      lambda = lambda,
      loglik = q[2],
      all.loglik = rev(q[-1])[-1],
      iter = it,
      size = size
    ),
    class = "gaussian_mixture"
  )
}

##############################################
### I don't want to export these functions ###
##############################################

# @title Density of a multivariate normal distribution
#
# @description Compute the density of a multivariate normal distribution
#
# @param x is a matrix where each row is a data point
# @param mu is a vector
# @param sigma is a square matrix with sides as big as x has columns
#
# @return a vector with the density of each data point
#
# @author Eduardo Ruiz Sabajanes, \email{eduardoruizsabajanes@@gmail.com}
dmnorm <- function(x, mu, sigma) {
  k <- ncol(sigma)

  x  <- as.matrix(x)
  diff <- t(t(x) - mu)

  num <- exp(-1 / 2 * diag(diff %*% solve(sigma) %*% t(diff)))
  den <- sqrt(((2 * pi)^k) * det(sigma))
  num / den
}

# Finite versions of sum, rowSums, colSums and cov.wt i.e. versions that
# ignore NA, NaN and Inf values
sum.finite <- function(x) {
  sum(x[is.finite(x)])
}

rowSums.finite <- function(x) {
  rowSums(x[, is.finite(colSums(x)), drop = FALSE])
}

colSums.finite <- function(x) {
  colSums(x[is.finite(rowSums(x)), , drop = FALSE])
}

cov.wt.finite <- function(x, wt, center) {
  fnt <- is.finite(wt)
  cov.wt(x[fnt, , drop = FALSE], wt = wt[fnt], center = center)
}
