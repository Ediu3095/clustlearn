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
distance <- function(c1, c2 = NULL, ..., simplify = TRUE) {
  d <- if (is.null(c2) && simplify) {
    dist(c1, ...)
  } else if (is.null(c2)) {
    as.matrix(dist(c1, ...))
  } else {
    rows <- seq.int(from = nrow(c1) + 1, length.out = nrow(c2))
    cols <- seq.int(length.out = nrow(c1))
    as.matrix(dist(rbind(c1, c2), ...))[rows, cols]
  }
  if (!is.null(c2) || !simplify) {
    colnames(d) <- rownames(c1)
    rownames(d) <- rownames(c2)
  }
  d
}
