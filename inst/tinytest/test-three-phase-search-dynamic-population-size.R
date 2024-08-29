library(anticlustPackage)
library(Rcpp)

Rcpp::evalCpp(2 + 2)
Rcpp::compileAttributes()

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

ergebnis <- anticlustPackage::three_phase_search_anticlustering(dat, K, N)

result <- .Call("_anticlustPackage_three_phase_search_dynamic_population_size", dat, N, K, upper_bound, lower_bound, popSize, time_limit, theta_max, theta_min, beta_min, LMAX)

.Call(`_anticlustPackage_anticlust_three_phase_search_dynamic_population_size`, matrix, N, K, upper_bound, lower_bound, popSize, time_limit, theta_max, theta_min, beta_min, LMAX)


library(devtools)

devtools::document()  # To update NAMESPACE and documentation
devtools::load_all()  # To load the package and make functions available in the current session

devtools::document()  # Update documentation and NAMESPACE
devtools::build()     # Build the package
devtools::install()   # Install the package

