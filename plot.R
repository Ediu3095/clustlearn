png(
  file = '...',
  width = 500,
  height = 500
)
par(mar = c(1, 1, 1, 1))
plot(..., hang = -1, main = '', sub = '', xlab = '', ylab = '', asp = 1, frame.plot = TRUE, axes = FALSE, pch = 20)
dev.off()

###############################################################################################

# VONOROI DIAGRAM

library(deldir)

x <- matrix(rnorm(1000, sd = 0.2), ncol = 2)
colnames(x) <- c("x", "y")
cl <- kmeans(x, matrix(c(-.2, -.2, .2, .2, -.2, .2, -.2, .2), ncol = 2), max_iterations = 0)

tesselation <- deldir(cl$centers[, 1], cl$centers[, 2], rw = c(-1, 1, -1, 1))
tiles <- tile.list(tesselation)

png(
  file = 'D:\\Pictures\\TFG\\stateoftheart\\03-partitional-clustering-example.png',
  width = 500,
  height = 500
)
par(mar = c(1, 1, 1, 1))
plot(x, col = cl$cluster, xlim = c(-.7, .7), ylim = c(-.7, .7), main = '', sub = '', xlab = '', ylab = '', asp = 1, frame.plot = TRUE, axes = FALSE, pch = 20)
plot(tiles, main = '', sub = '', xlab = '', ylab = '', asp = 1, frame.plot = TRUE, axes = FALSE, pch = 19, add = TRUE, showpoints = FALSE, border = "#ffffff00", fillcol = c("#00000019", "#DF536B19", "#61D04F19", "#2297E619"))
dev.off()

###############################################################################################

png(
  file = 'D:\\Pictures\\TFG\\development\\db6.png',
  width = 500,
  height = 500
)
par(mar = c(1, 1, 1, 1))
plot(
  read.csv('C:\\Users\\3arci\\Desktop\\papers\\databases\\db6.csv'),
  main = '',
  sub = '',
  xlab = '',
  ylab = '',
  asp = 1,
  frame.plot = TRUE,
  axes = FALSE,
  pch = 20
)
dev.off()

###############################################################################################

# NORMAL DISTRIBUTION
dmnorm <- function(x, mu, sigma) {
  k <- ncol(sigma)

  x  <- as.matrix(x)
  diff <- t(t(x) - mu)

  num <- exp(-1 / 2 * diag(diff %*% solve(sigma) %*% t(diff)))
  den <- sqrt(((2 * pi)^k) * det(sigma))
  num / den
}

data <- read.csv('C:\\Users\\3arci\\Desktop\\papers\\databases\\db4.csv')

k <- 3
members <- kmeans(data, k, 100)
mu <- t(members$centers)
sigma <- sapply(
  seq_len(k),
  function(i) cov(data[members$cluster == i, , drop = FALSE]),
  simplify = "array"
)
lambda <- members$size / nrow(data)

x <- seq(-2.5, 2.5, .05)
y <- seq(-2.5, 2.5, .05)

par(mar = c(1, 1, 1, 1))
plot(data, col = members$cluster, main = '', sub = '', xlab = '', ylab = '', asp = 1, frame.plot = TRUE, axes = FALSE, pch = 20)
for (i in seq_len(k)) {
  m <- mu[, i]
  s <- sigma[, , i]
  f <- function(x, y) dmnorm(cbind(x, y), m, s)
  z <- outer(x, y, f)
  contour(x, y, z, col = i, add = TRUE)
}

members <- gaussian_mixture(data, k, 100)
mu <- members$mu
sigma <- members$sigma
lambda <- members$lambda

par(mar = c(1, 1, 1, 1))
plot(data, col = members$cluster, main = '', sub = '', xlab = '', ylab = '', asp = 1, frame.plot = TRUE, axes = FALSE, pch = 20)
for (i in seq_len(k)) {
  m <- mu[, i]
  s <- sigma[, , i]
  f <- function(x, y) dmnorm(cbind(x, y), m, s)
  z <- outer(x, y, f)
  contour(x, y, z, col = i, add = TRUE)
}

plot(data, col = (dbscan(data, .25, 5) + 3) %% 4 + 1, main = '', sub = '', xlab = '', ylab = '', asp = 1, frame.plot = TRUE, axes = FALSE, pch = 20)

###############################################################################################
