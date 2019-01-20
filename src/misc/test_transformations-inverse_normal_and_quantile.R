
### Min/max scaling
x <- runif(n=100, min = -10, max = 100)
x
x.mm <- (x - min(x)) / (max(x) - min(x))


### Inverse Gaussian transformation
k <- 3/8
n <- 1000
x.raw <- runif(n, min = -10, max = 100)
qplot(x.raw, geom="density")
x.rank <- rank(x.raw)
x.calc <- (x.rank-k)/(n-2*k+1)
x.int <- qnorm(x.calc)
qplot(x.int, geom="density")

### Inverse Uniform transformation
n <- 1000
x.raw <- runif(n, min = -10, max = 100)
x.raw <- rnorm(n)
x.rank <- rank(x.raw)
x.calc <- x.rank/n
x.iut <- qunif(x.calc)
qplot(x.iut, geom="density")
df.x <- tibble(x.raw, x.rank, x.calc, x.iut)
