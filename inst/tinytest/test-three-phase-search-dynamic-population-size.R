library(anticlust)
library(tinytest)

convert_to_group_pattern <- function(vec) {
  # Create a mapping of unique values to group identifiers
  unique_values <- unique(vec)
  group_map <- setNames(seq_along(unique_values), unique_values)
  # Convert the original vector to the group identifiers
  as.numeric(factor(vec, levels = unique_values))
}

set.seed(123)

N <- 12
M <- 5
K <- 2
dat <- matrix(rnorm(N * M), ncol = M)
distances <- dist(dat)
distances2 <- anticlust:::convert_to_distances(dat) 

results1 <- anticlust:::three_phase_search_anticlustering(dat, K, N)
result_cluster1 <- results1$result
diversity1 <- diversity_objective(distances, result_cluster1)

result_cluster2 <- anticlustering(distances, K=K, method="local-maximum", repetitions = 10)
diversity2 <- diversity_objective(distances,  result_cluster2)

result_cluster3 <- optimal_anticlustering(distances, objective = "diversity", K=K, solver = "lpSolve")
diversity3 <- diversity_objective(distances,  result_cluster3)

expect_equal(diversity1, diversity2)
expect_equal(diversity2, diversity3)
expect_true(all(convert_to_group_pattern(result_cluster1) == convert_to_group_pattern(result_cluster2)))
expect_true(all(result_cluster3 == result_cluster2))

N <- 12
M <- 2
K <- 3
dat <- matrix(rnorm(N * M), ncol = M)
distances <- dist(dat)
distances2 <- anticlust:::convert_to_distances(dat) 

ergebnis <- anticlust:::three_phase_search_anticlustering(dat, K, N)
diversity_objective(distances, ergebnis$result)

ergebnis2 <- anticlustering(distances, K=K, method="local-maximum", repetitions = 10)
diversity_objective(distances, ergebnis2)

ergebnis3 <- optimal_anticlustering(distances, objective = "diversity", K=K, solver = "lpSolve")
diversity_objective(distances, ergebnis3)

print(ergebnis$result)
print(ergebnis2)
print(ergebnis3)

plot_clusters(dat, clusters = ergebnis$result,within_connection = TRUE, show_axes = TRUE)
plot_clusters(dat, clusters = ergebnis2, within_connection = TRUE,show_axes = TRUE)
plot_clusters(dat, clusters = ergebnis3, within_connection = TRUE, show_axes = TRUE)

