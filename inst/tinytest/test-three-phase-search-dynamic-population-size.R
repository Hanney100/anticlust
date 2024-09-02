library(anticlust)

set.seed(123)

N <- 12
M <- 5
K <- 2
dat <- matrix(rnorm(N * M), ncol = M)
upper_bound <- 2
lower_bound <- 2
theta_max  <- 1.2
theta_min  <- 0.1
beta_min  <- 2
LMAX <- 3
time_limit <- 3
popSize <- 15

ergebnis <- anticlust:::three_phase_search_anticlustering(dat, K, N)
