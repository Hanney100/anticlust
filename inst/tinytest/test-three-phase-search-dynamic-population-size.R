library(anticlustPackage)
library(Rcpp)

Rcpp::evalCpp(2 + 2)

set.seed(123)

N <- 12
M <- 5
K <- 2
dat <- matrix(rnorm(N * M), ncol = M)

result <- anticlustPackage:::three_phase_search_anticlustering(dat, K, N)

anticlustPackage:::rcpp_hello_world()
