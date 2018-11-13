#' Binary Data Simulation
#'
#' This function simulates multivariate binary data from a latent factor model.
#' @param n - number of sites.
#' @param p - number of species
#' @param mu.alpha - mean for alpha parameter: default is 0
#' @param sd.alpha - sd for alpha parameter: default is 1
#' @param mu.beta - mean for beta parameter: default is 0
#' @param sd.beta - sd for beta parameter: default is 1
#' @return mat - matrix representation of presence / absence
#' @return vec - vectorized representation of presence / absence
#' @return alpha.true - true values of alpha
#' @return beta.true - true values of alpha
#' @return theta.true - true values of alpha
#' @return z.true - true values of alpha
#' @export

simProbitY <- function(n, p, mu.alpha = 0, sd.alpha = 1, mu.beta = 1, sd.beta = 1 ){
  Y <- matrix(0,n,p)
  z.true <- matrix(stats::rnorm(2 * n),nrow = n, ncol=2)
  theta.true <- matrix(stats::rnorm(2 * p), nrow=p, ncol=2)
  theta.true[1,2] <- 0
  while(theta.true[1,1] < 0){
    theta.true[1,1] <- stats::rnorm(1)
  }
  while(theta.true[2,2] < 0){
    theta.true[2,2] <- stats::rnorm(1)
  }

  alpha.true <- stats::rnorm(n, mu.alpha, sd.alpha)
  beta.true <- stats::rnorm(p, mu.beta, sd.beta)

  for (i in 1:n){
    probs <- stats::pnorm(alpha.true[i] + beta.true + theta.true %*% z.true[i,])
    Y[i,] <- stats::rbinom(p,1,probs)
  }
  Y.vec <- as.numeric(Y)
  return(list(mat = Y, vec = Y.vec, alpha.true = alpha.true, beta.true = beta.true, theta.true = theta.true, z.true = z.true))
}
