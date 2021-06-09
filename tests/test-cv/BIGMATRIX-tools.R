
library(biganalytics)
install.packages('bigalgebra')

library(bigalgebra)

X <- matrix(rnorm(4 * 3), nrow = 4)
Xbig <- as.big.matrix(X)

beta <- matrix(rnorm(3 * 1), nrow = 3)
as.matrix(Xbig %*% beta)
X %*% beta

biganalytics::apply(Xbig, 2, sum)
colSums(X)

row.idx <- c(2, 4)
str(Xbig[1:2,,drop=F])
