

library("anticlustPackage")

# constructing selection matrix from clusterin vector and reversing works
for (M in 1:4) {
  for (K in 2:5) {
    for (i in 3:8) {
      N <- K * i
      clusters <- sample(rep(1:K, N/K))
      selected <- anticlustPackage:::selection_matrix_from_clusters(clusters)
      clusters2 <- anticlustPackage:::clusters_from_selection_matrix(selected)
      expect_true(all(anticlustPackage:::order_cluster_vector(clusters) == clusters2))
    }
  }
}
