#' Count Data Simulation
#'
#' This function simulates multivariate Poisson count data from a latent factor model.
#' @param n - number of sites.
#' @param p - number of species
#' @param mu.alpha - mean for alpha parameter: default is 0
#' @param sd.alpha - sd for alpha parameter: default is 1
#' @param mu.beta - mean for beta parameter: default is 0
#' @param sd.beta - sd for beta parameter: default is 1
#' @return mat - matrix representation of abundance
#' @return vec - vectorized representation of abundance
#' @return alpha.true - true values of alpha
#' @return beta.true - true values of alpha
#' @return theta.true - true values of alpha
#' @return z.true - true values of alpha
#' @export

simPoissonY <- function(n, p, mu.alpha = 0, sd.alpha = 1 , mu.beta = 0, sd.beta = 1){
  alpha.true <- stats::rnorm(n, mu.alpha, sd.alpha)
  beta.true <- stats::rnorm(p, mu.alpha, sd.alpha)
  z.true <- matrix(stats::rnorm(n*2),n,2)
  theta.true <- matrix(stats::rnorm(p*2), nrow=p, ncol=2)
  theta.true[1,2] <- 0
  while(theta.true[1,1] < 0){
    theta.true[1,1] <- stats::rnorm(1)
  }
  while(theta.true[2,1] < 0){
    theta.true[2,1] <- stats::rnorm(1)
  }

  Z.theta <- z.true %*% t(theta.true)
  mu.true <- rep(alpha.true, each = p) + rep(beta.true, n) + as.numeric(t(Z.theta))
  mu.matrix <- matrix(mu.true, n, p, byrow=T)
  Y <- stats::rpois(n*p,exp(mu.true))
  Y.matrix <- matrix(Y, n, p , byrow=T)
  return(list(mat = Y.matrix, vec = Y, alpha.true = alpha.true, beta.true = beta.true, theta.true = theta.true, z.true = z.true))
}
