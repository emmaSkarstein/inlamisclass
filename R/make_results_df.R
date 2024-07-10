#' Summarize results from inla_is or inla_mcmc
#'
#' @param mcmc_results results from inla_is or inla_mcmc
#' @param niter number of iterations for the MCMC
#' @param nburnin number of iterations for the burn-in
#'
#' @return A dataframe with summary statistics.
#' @export
#'
make_results_df <- function(mcmc_results, niter,
                                 nburnin = length(mcmc_results)-niter){
  #minimum of all marginal log likelihoods:
  max_LL <- max(do.call(rbind, lapply(mcmc_results, `[[`, "mlik")))

  # all marginal lls, but subtract the minimum first (for numerical stability):
  LL <- do.call(rbind, lapply(mcmc_results, `[[`, "mlik")) - max_LL

  trace <- (nburnin+1):niter
  # weights
  WW <- exp(LL[trace]) / sum(exp(LL[trace]))

  betas <-  do.call(rbind, lapply(mcmc_results, `[[`, "betax"))[trace]
  betas0.025 <-  do.call(rbind, lapply(mcmc_results, `[[`, "betax0.025"))[trace]
  betas0.975 <-  do.call(rbind, lapply(mcmc_results, `[[`, "betax0.975"))[trace]

  betax_mean <- sum(betas*WW)
  betax_0.025 <- sum(betas0.025*WW)
  betax_0.975 <- sum(betas0.975*WW)

  dd_results_beta <- data.frame(mean = betax_mean,
                                "0.025quant" = betax_0.025,
                                "0.975quant" = betax_0.975)
  return(dd_results_beta)
}
