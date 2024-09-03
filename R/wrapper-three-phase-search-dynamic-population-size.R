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
#' @param K How many anticlusters should be created. 
#' @param popSize beta_max from paper
#' @param beta_min
#' @param lower_bound
#' @param upper_bound
#' @param theta_min
#' @param theta_max
#' @param time_limit
#' @param return
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
#' @note
#' 
#' @references
#' 
#' Xiao Yang et al. “A three-phase search approach with dynamic population size for solving # nolint
#' the maximally diverse grouping problem”. In: European Journal of Operational Research
#' 302.3 (2022). [SOURCE-CODE: https://raw.githubusercontent.com/toyamaailab/toyamaailab.github.io/main/resource/TPSDP_Code.zip],
#'  pp. 925–953. ISSN: 0377-2217. DOI: https://doi.org/10.1016/j.ejor.2022.02.003. 
#' 
three_phase_search_anticlustering <- function(x, K, N,
    upper_bound  = NULL, lower_bound  = NULL, popSize = 15, time_limit  = NULL, theta_max = NULL, theta_min = NULL, beta_min = NULL, LMAX=3) {

    #input_validation_threephase_search(matrix, K)
    distances <- convert_to_distances(x) 
    cat("Current distances:", distances, "\n")
  
    cat("Current N:", N, "\n")
    cat("Current K:", K, "\n")
    
    if (is.null(x)) {
      cat("x is NULL\n")
    } else {
      cat("Current x:", x, "\n")
    }
    
 
    if (is.null(lower_bound)) {
       lower_bound <-  round(N/K)
    } 
    if (is.null(upper_bound)) {
       upper_bound <-  round(N/K)
    } 
     cat("Current Upper Bound:", upper_bound, "\n")
     cat("Current Lower Bound:", lower_bound, "\n")
    
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

     # create empty matrix for results to use in C
     result_vector = numeric(N)
     cat("Current result_vector:", result_vector, "\n")
     
     results <- .C("three_phase_search_dynamic_population_size",
                  distances = as.double(distances),
                  N_in = as.integer(N),
                  K_in = as.integer(K),
                  upper_bound = as.integer(upper_bound),
                  lower_bound = as.integer(lower_bound),
                  Beta_max = as.integer(popSize),
                  time_limit = as.integer(time_limit),
                  Theta_max = as.double(theta_max),
                  Theta_min = as.double(theta_min),
                  Beta_min = as.integer(beta_min),
                  Lmax = as.integer(LMAX),
                  result = as.integer(result_vector),
                  cost = as.double(0.0),
                  mem_error = as.integer(0),
                  PACKAGE = "anticlust"
     )
     print(results[["mem_error"]])
     print(results[["cost"]])
     cat("Cualculated result_vector:", result_vector, "\n")
     
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