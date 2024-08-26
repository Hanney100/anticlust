#' Three phase search with dynamic population size heuristic
#'
#' This function implements the three phase search algorithm TPSPD for
#' anticlustering by Yang et al. (2022; <doi.org/10.1016/j.ejor.2022.02.003>).
#' The description of their algorithm is
#' given in Section 2 of their paper (in particular, see the
#' Pseudocode in Algorithm 1).
#' 
#' 
#' @param x The data input. Currently just a vector.
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
library(Rcpp)

three_phase_search_anticlustering <- function(x, K, popSize = 15, LMAX=3,
    upper_bound  = NULL, lower_bound  = NULL, time_limit  = NULL, theta_max = 15, theta_min = NULL, beta_min = NULL) {

    N <- length(x)
    if (is.null(upper_bound)) {
        upper_bound = K
    } 
    if (is.null(lower_bound)) {
       lower_bound = K
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
       	if (N <= 120) { Time_limit <- 3 }
        else if (N <= 240) { Time_limit <- 20 }
        else if (N <= 480) { Time_limit <- 120 }
        else if (N <= 960) { Time_limit <- 600 }
        else if (N <= 2000) { Time_limit <- 1200 }
        else if (N <= 3000) { Time_limit <- 3000 }
        else { Time_limit  <- 5000 }
    } 


  return(three_phase_search_dynamic_population_size(D, N, K, upper_bound, lower_bound, popSize, time_limit, theta_max, theta_min, beta_min, LMAX))

}