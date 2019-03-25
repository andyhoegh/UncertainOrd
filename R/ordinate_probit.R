#' Probit Ordination
#'
#' Fits latent factor model to binary data with probit link
#' @param num.mcmc - number of mcmc iterations
#' @param Y.vec - an n X p vectorized representation of presence data
#' @return alpha.samples
#' @return beta.samples
#' @return theta.samples
#' @return z.samples
#' @export

ordinate_probit <- function(num.mcmc, Y.vec){
  # Initialize Storage
  alpha.samples <- matrix(0, num.mcmc, n, byrow=T)
  beta.samples <- matrix(0,num.mcmc, p, byrow=T)
  theta.samples <- array(0, dim = c(num.mcmc,p,2))
  z.samples <- array(0, dim = c(num.mcmc,n,2))

  latent.z.samples <- array(0, dim = c(num.mcmc,n,p))

  # Prior Values
  mu.beta <- mu.alpha <-  0 # prior mean
  V.beta <- V.alpha <- 1 # prior variance
  P.beta <- P.alpha <- solve(V.beta)
  V.theta <- 1
  P.theta <- diag(1/V.theta, nrow=2, ncol=2)

  # Precompute values
  S.beta <- solve(n + P.beta)
  S.alpha <- solve(p + P.alpha)

  # Set up Thresholds for Truncated Normal
  lower.vals <- rep(-Inf, n * p)
  upper.vals <- rep(Inf, n * p)
  lower.vals[Y.vec == 1] <- 0
  upper.vals[Y.vec == 0] <- 0

  for (iter in 2:num.mcmc){
    # Vectorized draw for latent z
    mean.z <- rep(alpha.samples[iter - 1, ],p) + rep(beta.samples[iter - 1,], each=n) + as.numeric(z.samples[iter-1,,] %*% t(theta.samples[iter -1,,]))
    latent.z.samples[iter, , ] <- truncnorm::rtruncnorm(n = n * p, a=lower.vals, b= upper.vals, mean = mean.z, sd = 1)

    # sample beta
    for (j in 1:p){
      m.beta <- S.beta %*% (sum(latent.z.samples[iter,,j] - alpha.samples[iter - 1, ] - z.samples[iter-1,,] %*% theta.samples[iter-1,j ,]) + P.beta * mu.beta)
      beta.samples[iter,j] <- stats::rnorm(1,m.beta, sd = sqrt(S.beta))
    }

    # sample alpha
    for (i in 1:n){
      m.alpha <- S.alpha %*% (sum(latent.z.samples[iter,i,] - beta.samples[iter,] -  theta.samples[iter - 1, ,] %*% z.samples[iter-1, i,]) + P.alpha * mu.alpha)
      alpha.samples[iter,i] <- stats::rnorm(1,m.alpha, sd = sqrt(S.alpha))
    }


    # sample theta, which is p x 2
    # for identifiability, the upper diagonal elements should be zero,
    # and diagonal elements are strictly positive
    for (j in 1:p){
      if (j == 1){
        S.theta <- solve(t(z.samples[iter-1,,1])  %*% z.samples[iter-1,,1] + solve(V.theta))
        m.theta <- S.theta * z.samples[iter-1,,1] %*% (latent.z.samples[iter,,j] - alpha.samples[iter,] - beta.samples[iter,j])
        theta.samples[iter,j,1] <- truncnorm::rtruncnorm(n=1, a=0, b=Inf, mean = m.theta, sd = sqrt(S.theta))
        theta.samples[iter,j,2] <- 0
      } else if (j == 2){
        S.theta <- solve(t(z.samples[iter-1,,]) %*% z.samples[iter-1,,] + P.theta)
        m.theta <- S.theta %*% t(z.samples[iter-1,,]) %*% (latent.z.samples[iter,,j] - alpha.samples[iter,] - beta.samples[iter,j])
        theta.samples[iter,j,] <- truncnorm::rtruncnorm(n=1, a = c(-Inf,0), b = c(Inf, Inf), mean = m.theta, sd = sqrt(diag(S.theta)))
      } else {
        S.theta <- solve(t(z.samples[iter-1,,]) %*% z.samples[iter-1,,] + P.theta)
        m.theta <- S.theta %*% t(z.samples[iter-1,,]) %*% (latent.z.samples[iter,,j] - alpha.samples[iter,] - beta.samples[iter,j])
        theta.samples[iter,j,] <- mnormt::rmnorm(n=1, mean = m.theta, varcov = S.theta)
      }
    }

    # sample z, which is n x 2
    S.z <- solve(t(theta.samples[iter,,]) %*%  theta.samples[iter,,] + diag(2))
    for (i in 1:n){
      m.z <- S.z %*% t(theta.samples[iter,,]) %*% (latent.z.samples[iter,i,] - alpha.samples[iter,i] - beta.samples[iter,])
      z.samples[iter,i,] <- mnormt::rmnorm(n=1, mean = m.z, varcov = S.z)
    }
  }
  return(list(alpha.samples = alpha.samples, beta.samples = beta.samples, theta.samples=theta.samples, z.samples = z.samples))
}
