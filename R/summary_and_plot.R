#' Summarize results from inla_is or inla_mcmc
#'
#' @param mcmc_results results from inla_is or inla_mcmc
#' @param niter number of iterations for the MCMC
#' @param nburnin number of iterations for the burn-in
#'
#' @return A dataframe with summary statistics.
#' @export
make_results_df <- function(mcmc_results, niter = length(mcmc_results),
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

  means <- calculate_summary_statistics(all.summary.fixed, "mean", WW_vec)
  lower_quant <- calculate_summary_statistics(all.summary.fixed, "0.025quant", WW_vec)
  upper_quant <- calculate_summary_statistics(all.summary.fixed, "0.975quant", WW_vec)

  summary_moi <- data.frame(means,
                           "0.025quant" = lower_quant$"0.025quant",
                           "0.975quant" = upper_quant$"0.975quant")
  colnames(summary_moi) <- c("variable", "mean", "0.025quant", "0.975quant")

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
calculate_summary_statistics <- function(summary_df, summary_stat, WW_vec){
  weighted_avg_summary <- dplyr::filter(summary_df, .data$summary_statistic == summary_stat) |>
    dplyr::mutate(WW = WW_vec) |>
    dplyr::group_by(.data$variable) |>
    dplyr::summarize(statistic = sum(.data$value*.data$WW))

  colnames(weighted_avg_summary) <- c("variable", summary_stat)
  return(weighted_avg_summary)
}

#' Plot inlamisclass models
#'
#' @inheritParams make_results_df
#'
#' @return A ggplot2 object.
#' @export
#' @importFrom ggplot2 ggplot aes vars
#' @importFrom rlang .data
plot_inlamisclass <- function(mcmc_results, niter = length(mcmc_results),
                              nburnin = length(mcmc_results)-niter){
  summary_df <- make_results_df(mcmc_results, niter = niter, nburnin = nburnin)
  ggplot2::ggplot(summary_df$moi, aes(y = .data$variable)) +
    ggplot2::geom_point(aes(x = .data$mean)) +
    ggplot2::geom_errorbarh(aes(xmin = .data$"0.025quant",
                                xmax = .data$"0.975quant"),
                            height=.2) +
    ggplot2::theme_bw()

}

#' Plot and compare inlamisclass models with other INLA models
#'
#' @inheritParams make_results_df
#' @param mcmc_results an inlamisclass model or list of inlamisclass models.
#' @param naive_mod a model that does not account for misclassification
#' @param plot_intercept should estimated intercept be plotted?
#' @param num_inlamisclass_models number of inlamisclass models if mcmc_results is a list
#'
#' @return A ggplot2 object.
#' @export
#' @importFrom ggplot2 ggplot aes vars
#' @importFrom rlang .data
plot_compare_inlamisclass <- function(mcmc_results, niter = NULL,
                              nburnin = 0, naive_mod, plot_intercept = FALSE, num_inlamisclass_models = 1){

  if(is.null(niter)){
    niter <- length(mcmc_results)
  }
  if(num_inlamisclass_models == 1){
    mcmc_results <- list(mcmc_results)
  }
  summary_list <- list()
  for(i in 1:num_inlamisclass_models){
    summary_list[[i]] <- make_results_df(mcmc_results[[i]], niter = niter, nburnin = nburnin)$moi
  }
  summary_inlamisclass <- dplyr::bind_rows(summary_list, .id = "submodel")
  summary_naive <- naive_mod$summary.fixed
  summary_naive$variable <- rownames(summary_naive)

  if(!plot_intercept){
    summary_inlamisclass <- dplyr::filter(summary_inlamisclass, .data$variable != "(Intercept)")
    summary_naive <- dplyr::filter(summary_naive, .data$variable != "(Intercept)")
  }

  if(num_inlamisclass_models == 1){summary_inlamisclass$submodel <- rep("", nrow(summary_inlamisclass))}
  summary_naive$submodel <- ""

  all_summary <- dplyr::bind_rows(inlamisclass = summary_inlamisclass, naive = summary_naive, .id = "model_type")

  all_summary$model <- paste0(all_summary$model_type, "[", all_summary$submodel, "]")

  all_summary$greek <- paste0("beta", "[", all_summary$variable, "]")

  ggplot2::ggplot(all_summary, aes(y = .data$model)) +
    ggplot2::geom_point(aes(x = .data$mean)) +
    ggplot2::geom_errorbarh(aes(xmin = .data$"0.025quant", xmax = .data$"0.975quant"), height=.2) +
    ggplot2::scale_y_discrete(limits = rev, labels = scales::parse_format()) +
    ggplot2::facet_wrap(vars(.data$greek), scales = "free_x", labeller = ggplot2::label_parsed) +
    ggplot2::theme_bw() +
    ggplot2::xlab("Posterior mean and 95% credible intervals") +
    ggplot2::theme(axis.title.y = ggplot2::element_blank())

}

