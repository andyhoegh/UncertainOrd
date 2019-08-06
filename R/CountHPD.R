#' Simulate HPD Width of Abundance Data
#'
#' function to simulate abundance data and calculated HPD width for latent variables
#' @param n - number of sites.
#' @param p - number of species
#' @param mu.alpha - mean for alpha parameter: default is 0
#' @param sd.alpha - sd for alpha parameter: default is 1
#' @param mu.beta - mean for beta parameter: default is 0
#' @param sd.beta - sd for beta parameter: default is 1
#' @return width - average HPD across all latent factors
#' @importFrom magrittr "%>%"
#' @examples
#' CountHPD(10,10)
#' @export

CountHPD <- function(n, p, mu.alpha = 0, sd.alpha = 1, mu.beta = 0, sd.beta = 1){
  sim.vals <- simPoissonY(n = n, p = p, mu.alpha = mu.alpha, sd.alpha = sd.alpha, mu.beta = mu.beta, sd.beta = sd.beta)
  count.mcmc <- boral::boral(sim.vals$mat, family = 'poisson', lv.control = list(num.lv = 2), row.eff = 'fixed', save.model =T, prior.control = list(type=c('normal','normal','normal','halfcauchy'), hypparams = c(1,1,1,20) ))
  set.seed(NULL) # as boral sets seed
  count.mcmc.sims <- count.mcmc$jags.model$BUGSoutput$sims.matrix
  count.mcmc.sims1  <- count.mcmc.sims %>% as.data.frame %>% dplyr::select(dplyr::starts_with('lvs[')) %>% dplyr::select(dplyr::ends_with('1]'))
  count.mcmc.sims2 <- count.mcmc.sims %>% as.data.frame %>% dplyr::select(dplyr::starts_with('lvs[')) %>% dplyr::select(dplyr::ends_with('2]'))
  count.mcmc.interval <- coda::HPDinterval(coda::as.mcmc(cbind(count.mcmc.sims1, count.mcmc.sims2)))
  count.mcmc.width <- round(mean(abs(count.mcmc.interval[,'upper'] - count.mcmc.interval[,'lower'])),2)
  return(width = count.mcmc.width )
}
