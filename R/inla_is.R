#' Fit model with misclassified covariate with INLA and importance sampling
#'
#' @inheritParams inla_is_mcmc
#' @param alpha values for coefficients of imputation model
#' @param ncores number of cores for the parallel computation.
#' @param conditional if the misclassification probabilities are conditional on some variable, then this should be the name of that variable.
#'
#'
#' @return A data frame with the resulting model.
#' @export
#'
inla_is <- function(formula_moi, formula_imp = NULL,
                    alpha, MC_matrix, data, niter, ncores = 4,
                    conditional = NULL, ...){
  r.out.naive <- INLA::inla(formula_moi, data = data, ...)

  models_list <- list()
  models_list[[1]] <- r.out.naive

  n <- nrow(data)

  vars <- identify_vars(formula_moi = formula_moi, formula_imp = formula_imp)

  mc.out <- parallel::mclapply(1:niter, function(i) {

    if(is.null(conditional)){
      sample_pi <- new_pi(alpha = alpha, z = data[, vars$imp_covs],
                          MC_matrix, w = data[, vars$error_var])
    }else{
      sample_pi <- new_pi_conditional(alpha = alpha, z = data[, vars$imp_covs],
                                      MC_matrix, w = data[, vars$error_var],
                                      conditioning_var = data[, conditional])
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
                     #control.inla = list(int.strategy = 'eb'),
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


#' Fit model with misclassified covariate with INLA and iterative IS
#'
#' @param formula_moi formula for the model of interest
#' @param formula_imp formula for the imputation model
#' @param alpha initial values for coefficients for imputation model
#' @param MC_matrix misclassification matrix
#' @param data data for the model
#' @param niter number of iterations
#' @param ... further arguments to be passed to inla().
#'
#' @return A data frame with the resulting model.
#' @export
#'
inla_is_iterative <- function(formula_moi, formula_imp = NULL,
                      alpha, MC_matrix,
                      data, niter = 1600,
                      ...){
  # So far, this assumes the exposure model just has one covariate.
  # And the formula_moi does not really do much either

  r.out.naive <- INLA::inla(formula_moi, data = data, ...)

  models_list <- list()
  models_list[[1]] <- r.out.naive

  n <- nrow(data)

  mc.out <- list()

  vars <- identify_vars(formula_moi = formula_moi, formula_imp = formula_imp)

  for (ii in 1:niter){
    #print(sprintf("%d", ii))

    sample_pi <- new_pi(alpha, z = data[, vars$imp_covs],
                        MC_matrix, w = data[, vars$error_var])

    # We have a model for x; use the sample_pi probabilities derived above to sample from it
    xstar <- stats::rbinom(n, 1, sample_pi)

    new_data <- cbind(data, xstar = xstar)
    new_formula_moi <- stats::reformulate(response = vars$response, termlabels = c("xstar", vars$error_free_covs))

    if(ii < niter){
      r.inla <- INLA::inla(new_formula_moi,
                           data = new_data,
                           control.mode = list(result = models_list[[ii]],
                                               restart = TRUE),
                           control.compute = list(return.marginals = FALSE),
                           #control.inla = list(int.strategy = 'eb'),
                           ...)
    }else if(ii == niter){
      r.inla <- INLA::inla(new_formula_moi,
                           data = new_data,
                           control.mode = list(result = models_list[[ii]],
                                               restart = TRUE),
                           ...)
    }

    models_list[[ii+1]] <- r.inla

    # sample new alpha for the exposure model of x. To this end, use xstar to fit the model and then extract alpha
    new_formula_imp <- stats::reformulate(response = "xstar", termlabels = vars$imp_covs)
    r.inla.x <- INLA::inla(new_formula_imp, data = new_data,
                           family = "binomial",
                           control.compute = list(config = TRUE,
                                                  return.marginals = FALSE))
                           #control.inla = list(int.strategy = 'eb'))

    # This contains all n predicted X-values, and the estimates for alpha.0 and alpha.z
    r.alpha <- INLA::inla.posterior.sample(n = 1, r.inla.x)
    alpha <- r.alpha[[1]]$latent[c(n+1, nrow(r.alpha[[1]]$latent))]
    names(alpha) <- c("alpha.0", paste0("alpha.", vars$imp_covs))

    summary.fixed <- r.inla$summary.fixed
    variables <- rownames(summary.fixed)

    # Change name of xstar to error variable
    error_var_index <- which(variables == "xstar")
    variables[error_var_index] <- vars$error_var

    summary.fixed$variable <- variables

    results <- list(mlik = r.inla$mlik[1,1],
                    summary.fixed = summary.fixed,
                    alpha = alpha)

    mc.out[[ii]] <- results
  }

  return(mc.out)
}
