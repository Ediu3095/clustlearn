#' @title Gaussian mixture model
#'
#' @description Perform Gaussian mixture model clustering on a data matrix.
#'
#' @param data a set of observations, presented as a matrix-like object where
#' every row is a new observation.
#' @param k the number of clusters to find.
#' @param max_iter the maximum number of iterations to perform.
#' @param details a Boolean determining whether intermediate logs explaining how
#' the algorithm works should be printed or not.
#' @param waiting a Boolean determining whether the intermediate logs should be
#' printed in chunks waiting for user input before printing the next or not.
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
#'  \item{Expectation}{Compute how much is each observation expected to belong
#'  to each component of the GMM.}
#'  \item{Maximization}{Recompute the GMM according to the expectations from
#'  the E-step in order to maximize them.}
#' }
#'
#' The algorithm stops when the changes in the expectations are sufficiently
#' small or when a maximum number of iterations is reached.
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
#' ### !! This algorithm is very slow, so we'll only test it on some datasets !!
#'
#' ### Helper functions
#' dmnorm <- function(x, mu, sigma) {
#'   k <- ncol(sigma)
#'
#'   x  <- as.matrix(x)
#'   diff <- t(t(x) - mu)
#'
#'   num <- exp(-1 / 2 * diag(diff %*% solve(sigma) %*% t(diff)))
#'   den <- sqrt(((2 * pi)^k) * det(sigma))
#'   num / den
#' }
#'
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
#'     f <- function(x, y) cl$lambda[i] * dmnorm(cbind(x, y), m, s)
#'     z <- outer(x, y, f)
#'     contour(x, y, z, col = i, add = TRUE)
#'   }
#' }
#'
#' ### Example 1
#' test(clustlearn::db1, 2)
#'
#' ### Example 2
#' # test(clustlearn::db2, 2)
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
#' # test(clustlearn::db6, 3)
#'
#' ### Example 7 (with explanations, no plots)
#' cl <- clustlearn::gaussian_mixture(
#'   clustlearn::db5[1:20, ],
#'   3,
#'   details = TRUE,
#'   waiting = FALSE
#' )
#'
#' @author Eduardo Ruiz Sabajanes, \email{eduardo.ruizs@@edu.uah.es}
#'
#' @importFrom stats cov
#' @export
gaussian_mixture <- function(
  data,
  k,
  max_iter = 10,
  details = FALSE,
  waiting = TRUE,
  ...
) {
  data <- as.matrix(data)

  if (details) {
    hline()
    console.log("EXPLANATION:")
    console.log("")
    console.log("The Gaussian Mixture Model with Expectation Maximization (GMM with EM) algorithm aims to model the data as a Gaussian Mixture Model i.e. the weighted sum of several Gaussian distributions, where each component i.e. each Gaussian distribution, represents a cluster.")
    console.log("")
    console.log("The Gaussian distributions are parameterized by their mean vector (mu), covariance matrix (sigma) and mixing proportion (lambda). Initially, the mean vector is set to the cluster centers obtained by performing a k-means clustering on the data, the covariance matrix is set to the covariance matrix of the data points belonging to each cluster and the mixing proportion is set to the proportion of data points belonging to each cluster. The algorithm then optimizes the GMM using the EM algorithm.")
    console.log("")
    console.log("The EM algorithm is an iterative algorithm that alternates between two steps:")
    console.log("")
    console.log("    1. Expectation step. Compute how much is each observation expected to belong to each component of the GMM.")
    console.log("    2. Maximization step. Recompute the GMM according to the expectations from the E-step in order to maximize them.")
    console.log("")
    console.log("The algorithm stops when the changes in the expectations are sufficiently small or when a maximum number of iterations is reached.")
    console.log("")

    if (waiting) {
      invisible(readline(prompt = "Press [enter] to continue"))
      console.log("")
    }
  }

  # Perform k-means to get initial values for mu, sigma and pi
  members <- kmeans(data, k, ...)
  mu <- members$centers
  sigma <- array(0, dim = c(k, ncol(data), ncol(data)))
  for (i in seq_len(k)) {
    sigma[i, , ] <- as.matrix(cov(data[members$cluster == i, , drop = FALSE]))
  }
  lambda <- members$size / nrow(data)

  if (details) {
    hline()
    console.log("INITIALIZATION:")
    console.log("")
    console.log("The GMM is initialized by calling kmeans. The initial components are:")
    for (i in seq_len(k)) {
      console.log(paste0("*** Component #", i, " ***\n"))
      cat("mu:\n")
      print(mu[i, ])
      cat("sigma:\n")
      print(as.matrix(sigma[i, , ]))
      cat("lambda:\n")
      print(lambda[i])
      console.log("")
    }
    console.log("These initial components are then optimized using the EM algorithm.")
    console.log("")

    if (waiting) {
      invisible(readline(prompt = "Press [enter] to continue"))
      console.log("")
    }
  }

  # EM algorithm
  # Starting values of expected value of the log likelihood
  q <- c(
    sum_finite(
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
  if (details) {
    hline()
    console.log("EM ALGORITHM:")
    console.log("")
    console.log("To measure how much the expectations change at each step we will use the log likelihood. The log likelihood is the sum of the logarithm of the probability of the data given the model. The higher the log likelihood, the better the model.")
    console.log("")
    console.log("The current log likelihood is:")
    cat("loglik:\n")
    print(q[1])
    console.log("")

    if (waiting) {
      invisible(readline(prompt = "Press [enter] to continue"))
      console.log("")
    }
  }
  while (abs(diff(q[1:2])) >= 1e-6 && it < max_iter) {
    # E step
    comp <- sapply(
      seq_len(k),
      function(i) lambda[i] * dmnorm(data, mu[i, ], as.matrix(sigma[i, , ]))
    )
    comp_sum <- rowSums_finite(comp)
    p <- comp / comp_sum

    if (details) {
      hline()
      console.log(paste0("STEP #", it, ":"))
      console.log("")
      console.log("E-STEP:")
      console.log("")
      console.log("The expectation of each observation to belong to each component of the GMM is the following:")
      cat("Expectation:\n")
      print(round(p, 3))
      console.log("")

      if (waiting) {
        invisible(readline(prompt = "Press [enter] to continue"))
        console.log("")
      }
    }

    # M step
    lambda <- sapply(
      seq_len(k),
      function(i) sum_finite(p[, i]) / nrow(data)
    )
    for (i in seq_len(k)) {
      mu[i, ] <- colSums_finite(p[, i] * data) / sum_finite(p[, i])
    }
    for (i in seq_len(k)) {
      tmp <- wtcov_finite(data, wt = p[, i], center = mu[i, ])$cov
      sigma[i, , ] <- as.matrix(tmp)
    }

    if (details) {
      console.log("M-STEP:")
      console.log("")
      console.log("The new components are:")
      for (i in seq_len(k)) {
        console.log(paste0("*** Component #", i, " ***\n"))
        cat("mu:\n")
        print(mu[i, ])
        cat("sigma:\n")
        print(as.matrix(sigma[i, , ]))
        cat("lambda:\n")
        print(lambda[i])
        console.log("")
      }

      if (waiting) {
        invisible(readline(prompt = "Press [enter] to continue"))
        console.log("")
      }
    }

    # Compute new expected value of the log likelihood
    q <- c(sum(log(comp_sum)), q)
    it <- it + 1

    if (details) {
      console.log("The new log likelihood is:")
      cat("loglik:\n")
      print(q[1])
      console.log("")

      if (waiting) {
        invisible(readline(prompt = "Press [enter] to continue"))
        console.log("")
      }
    }
  }

  comp <- sapply(
    seq_len(k),
    function(i) lambda[i] * dmnorm(data, mu[i, ], as.matrix(sigma[i, , ]))
  )
  cluster <- apply(comp, 1, which.max)
  size <- as.integer(table(factor(cluster, levels = seq_len(k))))

  if (details) {
    hline()
    console.log("FINAL RESULTS:")
    console.log("")
    if (it == max_iter)
      console.log("The algorithm stopped because the maximum number of iterations was reached.")
    else
      console.log("The algorithm stopped because the change in the log likelihood was smaller than 1e-6.")
    console.log("")
    console.log("With the current GMM every observation is assigned to the cluster it is most likely to belong to. The final clusters are:")
    cat("Cluster assignments:\n")
    print(cluster)
    console.log("")

    hline()
    console.log("")
  }

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
# @author Eduardo Ruiz Sabajanes, \email{eduardo.ruizs@@edu.uah.es}
dmnorm <- function(x, mu, sigma) {
  k <- ncol(sigma)

  x  <- as.matrix(x)
  diff <- t(t(x) - mu)

  num <- exp(-1 / 2 * diag(diff %*% solve(sigma) %*% t(diff)))
  den <- sqrt(((2 * pi)^k) * det(sigma))
  num / den
}

# Finite versions of sum, rowSums, colSums and cov.wt i.e. versions that
# replace NA, NaN and Inf values with 1e-300 (ignoring them would lead to
# numerical errors)
sum_finite <- function(x) {
  x[!is.finite(x)] <- 1e-300
  sum(x)
}

rowSums_finite <- function(x) {
  x[!is.finite(x)] <- 1e-300
  rowSums(x)
}

colSums_finite <- function(x) {
  x[!is.finite(x)] <- 1e-300
  colSums(x)
}

#' @importFrom stats cov.wt
wtcov_finite <- function(x, wt, center) {
  wt[!is.finite(wt)] <- 1e-300
  cov.wt(x, wt = wt, center = center)
}
