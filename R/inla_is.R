#' Fit model with misclassified covariate with INLA and importance sampling
#'
#' @inheritParams inla_mcmc
#' @param ncores number of cores for the parallel computation.
#' @param conditional Is the misclassification matrix conditional on other variables (so far only the moi response is possible, but this could easily be changed)
#'
#' @return A data frame with the resulting model.
#' @export
#'
inla_is <- function(formula_moi, formula_imp = NULL,
                    alpha, MC_matrix, data, niter, ncores = 4,
                    conditional = FALSE, ...){
  r.out.naive <- INLA::inla(formula_moi, data = data, ...)

  models_list <- list()
  models_list[[1]] <- r.out.naive

  n <- nrow(data)

  vars <- identify_vars(formula_moi = formula_moi, formula_imp = formula_imp)

  mc.out <- parallel::mclapply(1:niter, function(i) {

    if(!conditional){
      sample_pi <- new_pi(alpha, z = data[, vars$imp_covs],
                          MC_matrix, w = data[, vars$error_var])
    }else{
      sample_pi <- new_pi_conditional(alpha, z = data[, vars$imp_covs],
                                      MC_matrix, w = data[, vars$error_var])
    }

    # We have a model for x; use the sample_pi probabilities derived above to sample from it
    xstar <- stats::rbinom(n, 1, sample_pi)

    new_data <- cbind(data, xstar = xstar)
    new_formula_moi <- stats::reformulate(response = vars$response,
                                          termlabels = c("xstar", vars$error_free_covs))

    if(i < niter){
      r.inla <- INLA::inla(new_formula_moi,
                     data = new_data,
                     num.threads = 1,
                     control.mode = list(result = r.out.naive,
                                         restart = TRUE),
                     control.compute = list(return.marginals = FALSE),
                     control.inla = list(int.strategy = 'eb'),
                     ...
      )
    }else if(i == niter){
      r.inla <- INLA::inla(new_formula_moi,
                     data = new_data,
                     num.threads = 1,
                     control.mode = list(result = r.out.naive,
                                         restart = TRUE),
                     ...)
    }

    summary.fixed <- r.inla$summary.fixed

    # Change name of xstar to error variable
    variables <- rownames(summary.fixed)
    error_var_index <- which(variables == "xstar")
    variables[error_var_index] <- vars$error_var

    summary.fixed <- dplyr::mutate(summary.fixed, variable = variables)

    results <- list(mlik = r.inla$mlik[1,1],
                    summary.fixed = summary.fixed,
                    alpha = alpha)

    results

  }, mc.cores = ncores)
}

#' Fit model with misclassified covariate with INLA and importance sampling, where misclassification is conditional on response value
#'
#' @inheritParams inla_mcmc
#' @param ncores number of cores for the parallel computation.
#' @param MC_0 misclassification matrix for y = 0
#' @param MC_1 misclassification matrix for y = 1
#'
#' @return A data frame with the resulting model.
#' @export
#'
inla_is_conditional <- function(formula_moi, formula_imp = NULL,
                    alpha, MC_0, MC_1, data, niter, ncores = 4, ...){
  r.out.naive <- INLA::inla(formula_moi, data = data, ...)

  models_list <- list()
  models_list[[1]] <- r.out.naive

  n <- nrow(data)

  vars <- identify_vars(formula_moi = formula_moi, formula_imp = formula_imp)

  mc.out <- parallel::mclapply(1:niter, function(i) {

    sample_pi <- new_pi_conditional(alpha, z = data[, imp_covs],
                        MC_0 = MC_0, MC_1 = MC_1, w = data[, error_var])

    # We have a model for x; use the sample_pi probabilities derived above to sample from it
    xstar <- stats::rbinom(n, 1, sample_pi)

    new_data <- cbind(data, xstar = xstar)
    new_formula_moi <- stats::reformulate(response = response, termlabels = c("xstar", error_free_covs))

    r.inla <- INLA::inla(new_formula_moi,
                         data = new_data,
                         num.threads = 1,
                         control.mode = list(result = r.out.naive,
                                             restart = TRUE),
                         control.compute = list(return.marginals = FALSE),
                         control.inla = list(int.strategy = 'eb'),
                         ...)

    summary.fixed <- r.inla$summary.fixed

    # Change name of xstar to error variable
    variables <- rownames(summary.fixed)
    error_var_index <- which(variables == "xstar")
    variables[error_var_index] <- error_var

    summary.fixed <- dplyr::mutate(summary.fixed, variable = variables)

    results <- list(mlik = r.inla$mlik[1,1],
                    summary.fixed = summary.fixed,
                    alpha = alpha)

    results

  }, mc.cores = ncores)
}
