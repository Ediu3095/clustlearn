#' @title Expectation Maximization (EM) Algorithm
#'
#' @description Cluster analysis to find clusters of arbitrary shape and is
#' robust to outliers.
#'
#' @param data a set of observations, presented as a matrix-like object where
#' every row is a new observation.
#' @param k the number of clusters to find.
#' @param max_iter the maximum number of iterations to perform.
#' @param ... additional arguments passed to [kmeans()].
#'
#' @details
#'
#' @author Eduardo Ruiz Sabajanes, \email{eduardo.ruizs@@edu.uah.es}
#'
#' @return A [gaussian_mixture()] object.
#'
#' @examples
#'
#' @export
gaussian_mixture <- function(data, k, max_iter = 10, ...) {
  data <- as.matrix(data)

  # Perform k-means to get initial values for mu, sigma and pi
  members <- kmeans(data, k, ...)
  mu <- t(members$centers)
  sigma <- sapply(
    seq_len(k),
    function(i) as.matrix(cov(data[members$cluster == i, , drop = FALSE])),
    simplify = "array"
  )
  lambda <- members$size / nrow(data)

  # EM algorithm
  # Starting values of expected value of the log likelihood
  q <- c(
    sum.finite(
      sapply(
        seq_len(k),
        function(i) {
          log(lambda[i]) + log(dmnorm(data, mu[, i], sigma[, , i]))
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
      function(i) lambda[i] * dmnorm(data, mu[, i], sigma[, , i])
    )
    comp_sum <- rowSums.finite(comp)
    p <- comp / comp_sum

    # M step
    lambda <- sapply(
      seq_len(k),
      function(i) sum.finite(p[, i]) / nrow(data)
    )
    mu <- sapply(
      seq_len(k),
      function(i) colSums.finite(p[, i] * data) / sum(p[, i])
    )
    sigma <- sapply(
      seq_len(k),
      function(i) cov.wt.finite(data, wt = p[, i], center = mu[, i])$cov,
      simplify = "array"
    )

    # Compute new expected value of the log likelihood
    q <- c(sum(log(comp_sum)), q)
    it <- it + 1
  }

  comp <- sapply(
    seq_len(k),
    function(i) lambda[i] * dmnorm(data, mu[, i], sigma[, , i])
  )
  cluster <- apply(comp, 1, which.max)
  size <- as.integer(table(cluster))

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

# x is a matrix where each row is a data point
# mu is a vector
# sigma is a square matrix with sides as big as x has columns
dmnorm <- function(x, mu, sigma) {
  k <- ncol(sigma)

  x  <- as.matrix(x)
  diff <- t(t(x) - mu)

  num <- exp(-1 / 2 * diag(diff %*% solve(sigma) %*% t(diff)))
  den <- sqrt(((2 * pi)^k) * det(sigma))
  num / den
}

# Finite versions of sum, rowSums, colSums and cov.wt
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
