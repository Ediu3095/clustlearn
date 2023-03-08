#' @title Distance Matrix Computation
#'
#' @description This function computes and returns the distance matrix computed
#' by using the specified distance measure to compute the distances between the
#' rows of two data matrices.
#'
#' @param c1 a set of observations, presented as a matrix-like object where
#' every row is a new observation.
#' @param c2 a set of observations, presented as a matrix-like object where
#' every row is a new observation. If \code{c2} is not specified, it is
#' assumed to be equal to \code{c1}.
#' @param p the exponent of the Minkowski distance.
#' @param w the weight of every dimension in a weighted distance computation.
#' @param method the definition of distance. It must be one of
#' \code{"euclidean"}, \code{"manhattan"}, \code{"minkowski"}, \code{"octile"},
#' \code{"chebyshev"}, \code{"canberra"} or \code{"binary"}.
#' @param simplify a logical indicating whether results should be simplified if
#' possible. This is, whether the distance matrix should be cast to a [dist()]
#' object whenever possible.
#'
#' @details Available distance measures are (written for two vectors \eqn{x} and
#' \eqn{y}):
#'
#' \describe{
#'  \item{\code{binary}:}{(aka \emph{asymmetric binary}): The vectors are
#'  regarded as binary bits, so non-zero elements are 'on' and zero elements are
#'  'off'. The distance is the \emph{proportion} of bits in which only one is on
#'  amongst those in which at least one is on.}
#'
#'  \item{\code{canberra}:}{\eqn{\sum_i\left|x_i-y_i\right|/\left(\left|x_i
#'  \right|+\left|y_i\right|\right)}. Terms with zero denominator are omitted
#'  from the sum and treated as if the values were missing.}
#'
#'  \item{\code{chebyshev}:}{Maximum distance between two components of x and y
#'  (supremum norm), \eqn{\max_i\left(\left|x_i-y_i\right|\right)}.}
#'
#'  \item{\code{euclidean}:}{Usual distance between two vectors (2 norm aka
#'  \eqn{L_2}), \eqn{\sqrt{\sum_i\left(x_i-y_i\right)^2}}}
#'
#'  \item{\code{manhattan}:}{Absolute distance between the two vectors (1 norm
#'  aka \eqn{L_1}), \eqn{\sum_i\left|x_i-y_i\right|}}
#'
#'  \item{\code{minkowski}:}{The \eqn{p}th root of the sum of the \eqn{p}th
#'  powers of the absolute differences of the components (\eqn{p} norm aka
#'  \eqn{L_p}), \eqn{\left(\sum_i\left|x_i-y_i\right|^p\right)^\frac{1}{p}}}
#'
#'  \item{\code{octile}:}{The usual octile distance between \eqn{p(x,y)} and
#'  \eqn{p'(x',y')} is: \eqn{\sqrt{2} \times \min\left(\Delta{x},\Delta{y}
#'  \right) + \left|\Delta{x}-\Delta{y}\right|} where \eqn{\Delta{x} = \left|
#'  p_x-p'_x\right|} and \eqn{\Delta{y} = \left|p_y-p'_y\right|}. This function
#'  extends this formula to higher dimensional spaces.}
#' }
#'
#' Missing values are allowed, and are excluded from all computations. If some
#' columns are excluded in calculating the distance, the sum is scaled up
#' proportionally to the number of columns used. If all pairs are excluded when
#' calculating a particular distance, the value is \code{NaN}.
#'
#' @author Eduardo Ruiz Sabajanes, \email{eduardo.ruizs@@edu.uah.es}
#'
#' @return Whenever possible a [dist()] object. In any other case a [matrix()]
#' with columns corresponding to the objects in \code{c1} and rows to the
#' objects in \code{c2}.
#'
#' @examples
#' x <- matrix(rnorm(100), nrow = 5)
#' distance(x)
#'
#' ## example of binary and canberra distances.
#' x <- c(0, 0, 1, 1, 1, 1)
#' y <- c(1, 0, 1, 1, 0, 1)
#' distance(rbind(x, y), method = "binary")
#' ## answer 0.4 = 2/5
#' distance(rbind(x, y), method = "canberra")
#' ## answer 2 * (6/5)
#'
#' ## Examples involving "Inf" :
#' ## 1)
#' x[6] <- Inf
#' (m2 <- rbind(x, y))
#' distance(m2, method = "binary")   # warning, answer 0.5 = 2/4
#' ## These all give "Inf":
#' stopifnot(Inf == distance(m2, method =  "euclidean"),
#'           Inf == distance(m2, method =  "chebyshev"),
#'           Inf == distance(m2, method =  "manhattan"))
#'
#' @export
distance <- function(c1, c2 = c1, p = 3, w = NULL, method = "euclidean", simplify = TRUE) {
  # Load the adequate distance function according to method
  method <- switch(
    method,
    binary    = binary_distance,
    canberra  = canberra_distance,
    chebyshev = chebyshev_distance,
    euclidean = euclidean_distance,
    manhattan = manhattan_distance,
    minkowski = function(x, y, adj) minkowski_distance(x, y, p, adj),
    octile    = octile_distance
  )

  # Stop if the specified method is not supported
  if (is.null(method))
    stop("invalid distance method")

  d <- apply(
    c1,
    1,
    function(x) {
      apply(
        c2,
        1,
        compute_distance,
        x = x,
        w = w,
        method = method,
        simplify = TRUE
      )
    },
    simplify = TRUE
  )
  colnames(d) <- rownames(c1)
  rownames(d) <- rownames(c2)
  if (identical(c1, c2) && simplify)
    d <- stats::as.dist(d)
  d
}

compute_distance <- function(x, y, w = NULL, method) {
  # Convert inputs to numeric vectors
  x <- as.numeric(x)
  y <- as.numeric(y)

  # Ensure both vectors have the same dimensionality
  if (length(x) != length(y))
    stop("The x and y vectors must have the same dimensionality")

  # Identify the columns with which to make calculations
  z <- is.finite(x) | is.finite(y)
  if (identical(method, canberra_distance))
    z <- z | (x == 0 & y == 0)

  # Apply weights if necessary
  if (!is.null(w)) {
    # First make sure the dimensionality of the weights is alright
    w <- as.numeric(w)
    if (length(x) != length(w))
      stop("The x, y and w vectors must have the same dimensionality")

    # Then apply them
    z <- z & !(is.na(w) | is.nan(w))
    x[z] <- x[z] * w[z]
    y[z] <- y[z] * w[z]
  }

  # Filter both input vectors
  x <- x[z]
  y <- y[z]

  # Compute the adjustment factor
  adj <- length(z) / sum(z)

  # Compute the distance using the specified method
  if (sum(z) == 0) NaN else method(x, y, adj)
}

binary_distance <- function(x, y, adj) {
  # Compute the Binary distance
  sum(xor(x, y)) / sum(x | y) # * adj
}

canberra_distance <- function(x, y, adj) {
  # Compute the adjusted Canberra distance
  sum(abs(x - y) / (abs(x) + abs(y))) * adj
}

chebyshev_distance <- function(x, y, adj) {
  # Compute the Chebyshev distance
  max(abs(x - y))
}

euclidean_distance <- function(x, y, adj) {
  # Compute the Euclidean distance
  sqrt(sum((x - y) ^ 2) * adj)
}

manhattan_distance <- function(x, y, adj) {
  # Compute the Manhattan distance
  sum(abs(x - y) * adj)
}

minkowski_distance <- function(x, y, p, adj) {
  # Compute the Minkowski distance
  sum((abs(x - y) ^ p) * adj) ^ (1 / p)
}

octile_distance <- function(x, y, adj) {
  # Compute the absolute difference of each dimension and sort them
  dif <- -sort(-abs(x - y))

  # Compute the difference between each dimension and the previous one
  dif <- dif - c(dif[-1], 0)

  # Compute the Octile distance
  sum(dif * sqrt(seq_along(dif)) * adj)
}
