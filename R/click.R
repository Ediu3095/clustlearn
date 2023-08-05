click <- function(data) {
  # Preprocess data
  data <- scale(as.matrix(data))

  # Basic CLICK Algorithm
  s <- data %*% t(data)

  gaussian_mixture(as.vector(s), 2)
  mean_T <- 0
  std_T  <- 0
  mean_F <- 0
  std_F  <- 0
}
