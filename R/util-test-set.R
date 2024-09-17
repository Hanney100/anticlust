#' Generate Random Clusters of Data Points
#'
#' This function generates random clusters of data points with specified
#' parameters, ensuring that clusters are sufficiently spaced apart to avoid overlap.
#'
#' @param N_clusters Integer. The number of clusters to generate. Default is 7.
#' @param std Double. The standard deviation for generating cluster points,
#'     controlling the spread of the clusters. Default is 0.2.
#' @param N_data_per_cluster Integer. The number of data points to generate per cluster.
#'     Default is 6.
#' @param seed Integer. Random seed for reproducibility. Default is 123.
#'
#' @return A list containing the following elements:
#' \item{data}{A matrix of shape 2 x (N_clusters * N_data_per_cluster), where each column 
#'     represents the coordinates (X, Y) of a data point.}
#' \item{dist_matrix}{The distance matrix, computed as the dot product of the data matrix
#'     and its transpose.}
#' \item{labels}{A vector of length (N_clusters * N_data_per_cluster) containing the
#'     cluster labels for each data point.}
#'
#' @examples
#' result <- generate_test_clusters(N_clusters = 5, std = 0.2, N_data_per_cluster = 10)
#' dim(result$data)
#' dim(result$labels)
#'
#' @noRd
generate_test_clusters <- function(N_clusters = 7, std = 0.2, N_data_per_cluster = 6, seed = 123) {
  set.seed(seed)  # for reproducibility
  
  # Generate cluster means
  cluster_means <- matrix(rnorm(N_clusters * 2, mean = 0, sd = 20 * std), ncol = 2)
  
  # Ensure that clusters are not too close
  for (i in 1:(N_clusters - 1)) {
    for (j in (i + 1):N_clusters) {
      stopifnot(sum(abs(cluster_means[i, ] - cluster_means[j, ])) > 3 * std)
    }
  }
  
  # Generate data points around cluster means
  data <- do.call(rbind, lapply(1:N_clusters, function(i) {
    matrix(rnorm(N_data_per_cluster * 2, mean = cluster_means[i, ], sd = std), ncol = 2)
  }))
  
  dist_matrix <- matrix(0, nrow = N_clusters * N_data_per_cluster, ncol = N_clusters * N_data_per_cluster)
  # Calculate pairwise Euclidean distances
  for (i in 1:(N_clusters * N_data_per_cluster)) {
    for (j in 1:(N_clusters * N_data_per_cluster)) {
      if (i == j) {
        next  # Skip diagonal elements
      }
      dist_matrix[i, j] <- sqrt(sum((data[i, ] - data[j, ])^2))  # Euclidean distance
    }
  }
  
  # Generate labels for each data point
  labels <- rep(1:N_clusters, each = N_data_per_cluster)
  
  return(list(data = data, dist_matrix = dist_matrix, labels = labels))
}

#' Plot Clusters of Data Points
#'
#' This function creates a 2D scatter plot of data points colored by cluster.
#' It is intended to visualize the results of the `generate_test_clusters` function.
#'
#' @param data A matrix of shape 2 x (N_clusters * N_data_per_cluster), where each column 
#'     represents the coordinates (X, Y) of a data point.
#' @param labels A numeric vector of length (N_clusters * N_data_per_cluster) containing 
#'     the cluster labels for each data point.
#' @param N_clusters Integer. The number of clusters to plot. Default is 7.
#' @param N_data_per_cluster Integer. The number of data points per cluster. Default is 6.
#' @param std Double. The standard deviation for generating cluster points,
#'     controlling the spread of the clusters. Default is 0.2.
#' @param save_file Boolean. Saves the clustering image to a file.
#'
#' @return This function does not return any value. It creates a scatter plot.
#'
#' @examples
#' result <- generate_test_clusters()
#' plot_test_clusters(result$data, result$labels)
#'
plot_test_clusters <- function(data, labels, N_clusters = 7, N_data_per_cluster = 6, std = 0.2, save_file=FALSE) {
  # Define colors for up to 7 clusters
  colors <- rep(c("red", "blue", "green", "black", "magenta", "cyan", "yellow")[1:N_clusters], each = N_data_per_cluster)
  
 plot_title <- sprintf("N_data_per_cluster = %d, N_clusters = %d", N_data_per_cluster , N_clusters)

  if (save_file == TRUE) {
    file_name <- sprintf("cluster_plot_ M%d_N%d_std%f.png", N_data_per_cluster, N_clusters, std)
    png(file_name)
    plot(data[,1 ], data[, 2], col = colors, pch = 19, xlab = "X", ylab = "Y", main = plot_title)
    dev.off()
    message(sprintf("Plot saved as '%s'", file_name))
  } else {
    plot(data[, 1], data[, 2], col = colors, pch = 19, xlab = "X", ylab = "Y", main = plot_title)
  }
}