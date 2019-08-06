#' HPD Width Simulation
#'
#' This function simulates multivariate data from a latent factor model to estimate the average width of the HPD intervals
#' @param family - family of the distribution, "binomial" and "poisson" are currently supported
#' @param n - number of sites.
#' @param p - number of species
#' @param mu.alpha - mean for alpha parameter: default is 0
#' @param sd.alpha - sd for alpha parameter: default is 1
#' @param mu.beta - mean for beta parameter: default is 0
#' @param sd.beta - sd for beta parameter: default is 1
#' @param num.replicates - number of simulations: default is 1
#' @return width - average HPD across all latent factors
#' @examples HPDwidth('binomial',10,10)
#' @export

HPDwidth <- function(family, n, p, num.replicates = 1 , mu.alpha = 0, sd.alpha = 1, mu.beta = 1, sd.beta = 1){
  if (family == 'binomial'){
    width = replicate(num.replicates, BinaryHPD(n,p,mu.alpha,sd.alpha,mu.beta,sd.beta), simplify = T)
  } else if(family == 'poisson'){
    width = replicate(num.replicates, CountHPD(n,p,mu.alpha,sd.alpha,mu.beta,sd.beta), simplify = T)
  } else {
    stop('family should be binomial or poisson')
  }
  return(width = width)
}
