#' Fit model with misclassified covariate with INLA and importance sampling
#'
#' @inheritParams inla_mcmc
#' @param ncores number of cores for the parallel computation.
#'
#' @return A data frame with the resulting model.
#' @export
#'
inla_is <- function(formula_moi, formula_imp = NULL,
                    alpha, MC_matrix, data, niter, ncores = 4){
  r.out.naive <- INLA::inla(formula_moi, data = data)

  models_list <- list()
  models_list[[1]] <- r.out.naive

  n <- nrow(data)

  # Identify response variable
  response <- all.vars(formula_moi)[1]

  # Identify error variable
  error_var <- all.vars(formula_imp)[1]

  # Identify error free covariates
  error_free_covs <- labels(stats::terms(formula_imp))

  mc.out <- parallel::mclapply(1:niter, function(i) {

    sample_pi <- new_pi(alpha, z = data[, error_free_covs],
                        MC_matrix, w = data[, error_var])

    # We have a model for x; use the sample_pi probabilities derived above to sample from it
    xstar <- stats::rbinom(n, 1, sample_pi)

    new_data <- cbind(data, xstar = xstar)
    new_formula_moi <- stats::reformulate(response = response, termlabels = c("xstar", error_free_covs))

    if(i < niter){
      r.inla <- INLA::inla(new_formula_moi,
                     data = new_data,
                     num.threads = 1,
                     control.mode = list(result = r.out.naive,
                                         restart = TRUE),
                     control.compute = list(return.marginals = FALSE),
                     control.inla = list(int.strategy = 'eb')
      )
    }else if(i == niter){
      r.inla <- INLA::inla(new_formula_moi,
                     data = new_data,
                     num.threads = 1,
                     control.mode = list(result = r.out.naive,
                                         restart = TRUE))
    }

    summary.fixed <- r.inla$summary.fixed
    summary.fixed$variable <- rownames(summary.fixed)

    results <- list(mlik = r.inla$mlik[1,1],
                    summary.fixed = summary.fixed,
                    alpha = alpha)

    results

  }, mc.cores = ncores)
}
