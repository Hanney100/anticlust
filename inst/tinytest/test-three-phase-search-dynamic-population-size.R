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

result_cluster1 <- anticlust:::three_phase_search_anticlustering(distances, K, N)
result_cluster2 <- anticlustering(distances, K=K, method="local-maximum", repetitions = 10)
result_cluster3 <- optimal_anticlustering(distances, objective = "diversity", K=K, solver = "lpSolve")

diversity1 <- diversity_objective(distances, result_cluster1$result)
diversity2 <- diversity_objective(distances, result_cluster2)
diversity3 <- diversity_objective(distances, result_cluster3)

expect_equal(diversity1, diversity3)
expect_equal(diversity2, diversity3)


result1 <- convert_to_group_pattern(result_cluster1$result)
result2 <- convert_to_group_pattern(result_cluster2)
result3 <- convert_to_group_pattern(result_cluster3)

expect_true(all(result1 == result3))
expect_true(all(result2 == result3))

result_cluster1$result
result1
result2
result3

### Test more clusters ###

N <- 12
M <- 2
K <- 3

dat <- matrix(rnorm(N * M), ncol = M)
distances <- dist(dat)

result_cluster1 <- anticlust:::three_phase_search_anticlustering(distances, K, N)
result_cluster2 <- anticlustering(distances, K=K, method="local-maximum", repetitions = 10)
result_cluster3 <- optimal_anticlustering(distances, objective = "diversity", K=K, solver = "lpSolve")


diversity1 <- diversity_objective(distances, result_cluster1$result)
diversity2 <- diversity_objective(distances, result_cluster2)
diversity3 <- diversity_objective(distances, result_cluster3)
result_cluster1$score

expect_equal(diversity1, diversity3, info = "tdspd and lpsolve diversity are identical")
expect_equal(diversity2, diversity3, info = "local maximum and lspsolve diversity are identical")


result1 <- convert_to_group_pattern(result_cluster1$result)
result2 <- convert_to_group_pattern(result_cluster2)
result3 <- convert_to_group_pattern(result_cluster3)

expect_true(all(result1 == result3), info = "tdspd and lpsolve are identical")
expect_true(all(result2 == result3), info = "local maximum and lspsolve are identical")

result_cluster1$result
result1
result2
result_cluster2
result3
result_cluster3

