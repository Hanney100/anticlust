library(anticlust)

set.seed(123)

N <- 12
M <- 5
K <- 2
dat <- matrix(rnorm(N * M), ncol = M)
distances <- dist(dat)
distances2 <- anticlust:::convert_to_distances(dat) 

ergebnis <- anticlust:::three_phase_search_anticlustering(dat, K, N)
diversity_objective(distances, ergebnis)

ergebnis2 <- anticlustering(distances, K=K, method="local-maximum", repetitions = 10)
diversity_objective(distances, ergebnis3)

ergebnis3 <- optimal_anticlustering(distances, objective = "diversity", K=K, solver = "gurobi")

diversity_objective(distances, ergebnis3)
print(ergebnis)
