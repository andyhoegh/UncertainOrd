#' Simulate HPD Width of Presence Data
#'
#' function to simulate presence data and calculated HPD width for latent variables
#' @param n - number of sites.
#' @param p - number of species
#' @param mu.alpha - mean for alpha parameter: default is 0
#' @param sd.alpha - sd for alpha parameter: default is 1
#' @param mu.beta - mean for beta parameter: default is 0
#' @param sd.beta - sd for beta parameter: default is 1
#' @return width - average HPD across all latent factors
#' @importFrom magrittr "%>%"
#' @examples
#' BinaryHPD(10,10)
#' @export

BinaryHPD <- function(n, p, mu.alpha = 0, sd.alpha = 1, mu.beta = 0, sd.beta = 1){
  sim.vals <- simProbitY(n = n, p = p, mu.alpha = mu.alpha, sd.alpha = sd.alpha, mu.beta = mu.beta, sd.beta = sd.beta)
  pres.mcmc <- boral::boral(sim.vals$mat, family = 'binomial', lv.control = list(num.lv = 2), row.eff = 'fixed', save.model =T, prior.control = list(type=c('normal','normal','normal','halfcauchy'), hypparams = c(1,1,1,20) ))
  set.seed(NULL) # as boral automatically sets seed
  pres.mcmc.sims <- pres.mcmc$jags.model$BUGSoutput$sims.matrix
  pres.mcmc.sims1  <- pres.mcmc.sims %>% as.data.frame %>% dplyr::select(dplyr::starts_with('lvs[')) %>% dplyr::select(dplyr::ends_with('1]'))
  pres.mcmc.sims2 <- pres.mcmc.sims %>% as.data.frame %>% dplyr::select(dplyr::starts_with('lvs[')) %>% dplyr::select(dplyr::ends_with('2]'))
  pres.mcmc.interval <- coda::HPDinterval(coda::as.mcmc(cbind(pres.mcmc.sims1, pres.mcmc.sims2)))
  pres.mcmc.width <- round(mean(abs(pres.mcmc.interval[,'upper'] - pres.mcmc.interval[,'lower'])),2)
  return(width = pres.mcmc.width)
}
