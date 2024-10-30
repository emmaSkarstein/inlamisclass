#' Summarize results from inla_is_misclass
#'
#' @param is_results results from inla_is_misclass
#' @param niter number of iterations for the importance sampling
#'
#' @return A dataframe with summary statistics.
#' @export
make_results_df <- function(is_results, niter = length(is_results)){

  # List of all log likelihoods
  LL_vec <- do.call(rbind, lapply(is_results, `[[`, "mlik"))
  # minimum of all marginal log likelihoods:
  max_LL <- max(LL_vec)
  # all marginal lls, but subtract the minimum first (for numerical stability):
  LL <- LL_vec - max_LL

  # weights
  WW <- exp(LL) / sum(exp(LL))

  num_variables <- nrow(is_results[[1]]$summary.fixed)
  WW_vec <- rep(WW, 1, each = num_variables)

  all.summary.fixed <- dplyr::bind_rows(
    lapply(is_results,
           `[[`, "summary.fixed"), .id = "iteration") |>
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

  summary_imp <- colMeans(do.call(rbind, lapply(is_results, `[[`, "alpha")))

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
plot_inlamisclass <- function(is_results, niter = length(is_results)){
  summary_df <- make_results_df(is_results, niter = niter)
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
#' @param naive_mod a model that does not account for misclassification
#' @param correct_mod a model that uses the true values of the covariate with misclassification (only possible in simulation studies)
#' @param plot_intercept should estimated intercept be plotted?
#' @param num_inlamisclass_models number of inlamisclass models if is_results is a list
#'
#' @return A ggplot2 object.
#' @export
#' @importFrom ggplot2 ggplot aes vars
#' @importFrom rlang .data
plot_compare_inlamisclass <- function(is_results, naive_mod, correct_mod = NULL,
                                      niter = NULL,
                                      plot_intercept = FALSE,
                                      num_inlamisclass_models = 1){

  if(is.null(niter)){
    niter <- length(is_results)
  }
  if(num_inlamisclass_models == 1){
    is_results <- list(is_results)
  }
  # Summary of misclassification model(s)
  summary_list <- list()
  for(i in 1:num_inlamisclass_models){
    summary_list[[i]] <- make_results_df(is_results = is_results[[i]],
                                         niter = niter)$moi
  }
  summary_inlamisclass <- dplyr::bind_rows(summary_list, .id = "submodel")

  # Summary of naive model
  summary_naive <- naive_mod$summary.fixed
  summary_naive$variable <- rownames(summary_naive)

  # Summary of correct model
  if(!is.null(correct_mod)){
    summary_correct <- correct_mod$summary.fixed
    summary_correct$variable <- rownames(summary_correct)

    w_name <- setdiff(summary_naive$variable, summary_correct$variable)
    x_name <- setdiff(summary_correct$variable, summary_naive$variable)

    # Rename to correct coef name (so they are all beta_x)
    summary_inlamisclass[which(summary_inlamisclass$variable == w_name), "variable"] <- x_name
    summary_naive[which(summary_naive$variable == w_name), "variable"] <- x_name
  }

  if(num_inlamisclass_models == 1){summary_inlamisclass$submodel <- rep("", nrow(summary_inlamisclass))}
  summary_naive$submodel <- ""


  if(is.null(correct_mod)){
    all_summary <- dplyr::bind_rows(inlamisclass = summary_inlamisclass,
                                    naive = summary_naive, .id = "model_type")
  }else{
    all_summary <- dplyr::bind_rows(inlamisclass = summary_inlamisclass,
                                    naive = summary_naive,
                                    correct = summary_correct, .id = "model_type")
  }

  if(!plot_intercept){
    all_summary <- dplyr::filter(all_summary, .data$variable != "(Intercept)")
  }

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


#' Calculate sample misclassification matrices from simulated data
#'
#' @param data a simulated data set with a misclassified covariate
#'
#' @return A ggplot2 object giving the sample misclassification matrix.
#' @export
#' @importFrom ggplot2 ggplot aes vars
#' @importFrom rlang .data
data_characteristics <- function(data){
  pix0 <- sum(data$x==0)
  pix1 <- sum(data$x==1)
  est_MC_matrix <- table(x = data$x, w = data$w)/matrix(c(pix0, pix0, pix1, pix1), byrow = TRUE, nrow = 2)

  # Convert misclassification matrix to data frame
  matrix_df <- as.data.frame(as.table(round(est_MC_matrix, 3)))

  # create confusion matrix with ggplot2
  p <- ggplot2::ggplot(matrix_df, aes(x = .data$w, y = .data$x)) +
    ggplot2::geom_tile(fill = "white", color = "black") +
    ggplot2::theme_bw() +
    ggplot2::coord_equal() +
    ggplot2::scale_y_discrete(limits = rev) +
    ggplot2::geom_text(aes(label = .data$Freq), color = "black") +
    # following lines only increase text size (optional)
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank())

  return(p)
}
