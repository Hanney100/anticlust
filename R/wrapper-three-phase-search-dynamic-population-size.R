#' Three phase search with dynamic population size heuristic
#'
#' This function implements the three phase search algorithm TPSPD for
#' anticlustering by Yang et al. (2022; <doi.org/10.1016/j.ejor.2022.02.003>).
#' The description of their algorithm is
#' given in Section 2 of their paper (in particular, see the
#' Pseudocode in Algorithm 1).
#' 
#' 
#' @param matrix The data input. Currently just a vector.
#' @param K Number of anticlusters to be formed.
#' @param clusters A vector of length K that specifies the number of elements each cluster can contain. If this vector is not NULL, the lower and upper bounds will be disregarded.
#' @param beta_max The algorithm begins with a pool of random initial solutions of size beta_max. 
#'  Over time, the size of the solution pool decreases linearly until it reaches beta_min.
#' @param beta_min The minimum solution pool size the algorithm should reach before making a determination.
#' @param lower_bound Minimum number of elements in each anticluster. By default, anticlusters are of equal size, calculated as the total number of items divided by the number of clusters.
#' @param upper_bound Maximum number of elements in each anticluster. By default, anticlusters are of equal size, calculated as the total number of items divided by the number of clusters.
#' @param theta_max Parameter for the strength of undirected perturbation, which decreases linearly over time from theta_max to theta_min.
#' @param theta_min Parameter for the strength of undirected perturbation, which decreases linearly over time from theta_max to theta_min.
#' @param eta_max Parameter for 
#' @param alpha Parameter for weitghing the discrimitation of a slighlty worse local optiomal child solution.
#' @param time_limit Maximum execution time of the algorithm (in seconds).
#' @param return results contains everzthing, including vector result and its cost
#'     
#' @details
#' 
#' @return 
#' 
#' @author Hannah Hengelbrock \email{Hannah.Hengelbrock@@hhu.de}, 
#' Martin Papenberg \email{martin.papenberg@@hhu.de}
#' 
#' @export
#' 
#' @examples 
#' 
#' # Generate some random data
#' N <- 12
#' M <- 5
#' K <- 2
#' dat <- matrix(rnorm(N * M), ncol = M)
#' distances <- dist(dat)
#'
#' # Perform three hase serach algorithm
#' ergebnis <- anticlust:::three_phase_search_anticlustering(dat, K, N)
#'
#' # Compute objectives funtion
#' diversity_objective(distances, ergebnis$ergebnis)
#' 
#' # Compute comparision function
#' ergebnis2 <- anticlustering(distances, K=K, method="local-maximum", repetitions = 10)
#' diversity_objective(distances, ergebnis3)
#' 
#' # Compare both results
#' print(ergebnis$result)
#' print(ergebnis2)
#' 
#' @note
#' Important! for windows and linux there is a differnt definition of thie run time due to clock(),
#' so the time_limit acts differnetly depending on zour operating system
#' On windows it is the wall time, on linux the CPU time
#' 
#' @references
#' 
#' Xiao Yang et al. “A three-phase search approach with dynamic population size for solving # nolint
#' the maximally diverse grouping problem”. In: European Journal of Operational Research
#' 302.3 (2022). [SOURCE-CODE: https://raw.githubusercontent.com/toyamaailab/toyamaailab.github.io/main/resource/TPSDP_Code.zip],
#'  pp. 925–953. ISSN: 0377-2217. DOI: https://doi.org/10.1016/j.ejor.2022.02.003. 
#' 
three_phase_search_anticlustering <- function(x, K, N,
    max_iterations=100, clusters = NULL,
    upper_bound  = NULL, lower_bound  = NULL, beta_max = 15, 
    time_limit  = NULL, theta_max = NULL, theta_min = NULL, 
    beta_min = NULL, eta_max=3, alpha=0.05) {

    #input_validation_threephase_search(x, K)
    distances <- convert_to_distances(x) 
 
    if (is.null(lower_bound)) {
       lower_bound <- floor(N/K)
    } 
    if (is.null(upper_bound)) {
       upper_bound <- ceiling(N/K)
    } 
    
    if (N <= 400  & is.null(theta_max) & is.null(theta_min) & is.null(beta_min)) {
    	theta_max  <- 1.2
    	theta_min  <- 0.1
    	beta_min  <- 2
    } else if ( is.null(theta_max) & is.null(theta_min) & is.null(beta_min)) {
    	theta_max  <- 2.0
    	theta_min  <- 1.0
    	beta_min  <- 1
    }

    if (is.null(time_limit)) {
       	if (N <= 120) { time_limit <- 3 }
        else if (N <= 240) { time_limit <- 20 }
        else if (N <= 480) { time_limit <- 120 }
        else if (N <= 960) { time_limit <- 600 }
        else if (N <= 2000) { time_limit <- 1200 }
        else if (N <= 3000) { time_limit <- 3000 }
        else { time_limit  <- 5000 }
    } 

     # create result vector for results to use in C
    result_vector <- numeric(N)
    
    # check if K and clusters match 
    if (is.null(clusters)) {
      clusters <- rep(-1, K)
    } else if (length(clusters) != K) {
       stop("Number of giving clusters is not K.")
    }
     
     results <- .C("three_phase_search_dynamic_population_size",
                  distances = as.double(distances),
                  N_in = as.integer(N),
                  K_in = as.integer(K),
                  number_of_iterations = as.integer(max_iterations),
                  clusters = as.integer(clusters),
                  upper_bound = as.integer(upper_bound),
                  lower_bound = as.integer(lower_bound),
                  Beta_max = as.integer(beta_max),
                  time_limit = as.integer(time_limit),
                  Theta_max = as.double(theta_max),
                  Theta_min = as.double(theta_min),
                  Beta_min = as.integer(beta_min),
                  Eta_max = as.integer(eta_max),
                  Alpha = as.double(alpha),
                  result = as.integer(result_vector),
                  score = as.double(0.0),
                  mem_error = as.integer(0),
                  PACKAGE = "anticlust"
     )

     results[["mem_error"]]
     if (results[["mem_error"]] == 1) {
       stop("Could not allocate enough memory.")
     }

    return(results)
}

input_validation_threephase_search <- function(x, K) {
  
  return(input_validation_anticlustering(
    x, K, objective = "diversity", method = "brusco", 
    preclustering = FALSE, categories = NULL,
    repetitions = 1, standardize = FALSE
  ))
}