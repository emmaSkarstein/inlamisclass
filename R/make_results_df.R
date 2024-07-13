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
  # List of all log likelihoods
  LL_vec <- do.call(rbind, lapply(mcmc_results, `[[`, "mlik"))
  # minimum of all marginal log likelihoods:
  max_LL <- max(LL_vec)
  # all marginal lls, but subtract the minimum first (for numerical stability):
  LL <- LL_vec - max_LL

  trace <- (nburnin+1):niter

  # weights
  WW <- exp(LL[trace]) / sum(exp(LL[trace]))

  num_variables <- nrow(mcmc_results[[1]]$summary.fixed)
  WW_vec <- rep(WW, 1, each = num_variables)

  all.summary.fixed <- dplyr::bind_rows(
    lapply(mcmc_results, `[[`, "summary.fixed"), .id = "iteration") |>
    tidyr::pivot_longer(
      cols = c("mean", "sd", "0.025quant", "0.5quant", "0.975quant", "mode"),
      names_to = "summary_statistic")

  means <- calculate_summary_statistics(all.summary.fixed, "mean")
  lower_quant <- calculate_summary_statistics(all.summary.fixed, "0.025quant")
  upper_quant <- calculate_summary_statistics(all.summary.fixed, "0.975quant")

  summary_moi <- data.frame(means,
                           "0.025quant" = lower_quant$"0.025quant",
                           "0.975quant" = upper_quant$"0.975quant")

  summary_imp <- colMeans(do.call(rbind, lapply(mcmc_results, `[[`, "alpha")))

  return(list(moi = summary_moi, imp = summary_imp))
}

#' Calculate summary statistics
#'
#' @param summary_df data frame containing summary.fixed for each iteration (for all fixed effects)
#' @param summary_stat the name of the summary statistic of interest
#'
#' @return a data frame with the weighted average across all iterations for each fixed effect of the model
#' @keywords internal
#' @importFrom rlang .data
calculate_summary_statistics <- function(summary_df, summary_stat){
  weighted_avg_summary <- dplyr::filter(summary_df, .data$summary_statistic == summary_stat) |>
    dplyr::mutate(WW = .data$WW_vec) |>
    dplyr::group_by(.data$variable) |>
    dplyr::summarize(statistic = sum(.data$value*.data$WW))

  colnames(weighted_avg_summary) <- c("variable", summary_stat)
  return(weighted_avg_summary)
}
