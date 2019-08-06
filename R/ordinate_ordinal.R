#' Ordinal Ordination
#'
#' Fits latent factor model to ordinal data with cumulative probit sampling model
#' @param num.mcmc - number of mcmc iterations
#' @param Y - an n X p matrix of ordinal data
#' @param burn.in - number of burn in samples, default is 100
#' @return alpha.samples
#' @return beta.samples
#' @return theta.samples
#' @return z.samples
#' @examples
#' set.seed(08052019)
#' Y <- matrix(base::sample(1:5, 50, replace = TRUE), 10, 5)
#' ordinate_ordinal(1000, Y)
#' @export

ordinate_ordinal <- function(num.mcmc, Y, burn.in = 100){
  # Initialize Storage
  mcmc <- boral::boral(Y, family = 'ordinal', lv.control = list(num.lv = 2), row.eff = 'fixed', save.model =T, prior.control = list(type=c('normal','normal','normal','halfcauchy'), hypparams = c(1,1,1,20) ))
  set.seed(NULL) # as boral automatically sets seed
  mcmc.sims <- mcmc$jags.model$BUGSoutput$sims.matrix
  z1  <- mcmc.sims %>% as.data.frame %>% dplyr::select(dplyr::starts_with('lvs[')) %>% dplyr::select(dplyr::ends_with('1]'))
  z2 <- mcmc.sims %>% as.data.frame %>% dplyr::select(dplyr::starts_with('lvs[')) %>% dplyr::select(dplyr::ends_with('2]'))
  alpha <- mcmc.sims %>% as.data.frame %>% dplyr::select(dplyr::starts_with('row.coefs.ID1')) %>% dplyr::select(dplyr::ends_with('1]'))
  beta <- mcmc.sims %>% as.data.frame %>% dplyr::select(dplyr::starts_with('lv.coefs[')) %>% dplyr::select(dplyr::ends_with('1]'))
  theta1 <- mcmc.sims %>% as.data.frame %>% dplyr::select(dplyr::starts_with('lv.coefs[')) %>% dplyr::select(dplyr::ends_with('2]'))
  theta2 <- mcmc.sims %>% as.data.frame %>% dplyr::select(dplyr::starts_with('lv.coefs[')) %>% dplyr::select(dplyr::ends_with('3]'))

  return(list(alpha.samples = alpha, beta.samples = beta, theta.samples= abind::abind(theta1, theta2, along = 3), z.samples = abind::abind(z1, z2, along = 3)))
}
