#' Fit model with misclassified covariate with INLA and importance sampling
#'
#' @inheritParams inla_mcmc
#' @param ncores number of cores for the parallel computation.
#'
#' @return A data frame with the resulting model.
#' @export
#'
inla_is <- function(formula_moi, formula_imp = NULL,
                    alpha, MC_matrix, data, niter, ncores = 4, ...){
  r.out.naive <- INLA::inla(formula_moi, data = data, ...)

  models_list <- list()
  models_list[[1]] <- r.out.naive

  n <- nrow(data)

  # Identify response variable
  response <- all.vars(formula_moi)[1]

  # Identify error variable
  error_var <- all.vars(formula_imp)[1]

  # Covariates in imputation model
  imp_covs <- labels(stats::terms(formula_imp))

  # Identify error free covariates
  error_free_covs <- setdiff(c(response, error_var), all.vars(formula_moi))

  mc.out <- parallel::mclapply(1:niter, function(i) {

    sample_pi <- new_pi(alpha, z = data[, imp_covs],
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
    variables[error_var_index] <- error_var

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

  # Identify response variable
  response <- all.vars(formula_moi)[1]

  # Identify error variable
  error_var <- all.vars(formula_imp)[1]

  # Covariates in imputation model
  imp_covs <- labels(stats::terms(formula_imp))

  # Identify error free covariates
  error_free_covs <- setdiff(c(response, error_var), all.vars(formula_moi))

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

#' Calculation of new probabilities for x, given the alpha (and rest)
#'
#' @param alpha coefficients for the imputation model
#' @param z vector of covariates for imputation model, and that the misclassification matrix is conditioned on.
#' @param MC_1 misclassification matrix for z = 1
#' @param MC_0 misclassification matrix for z = 0
#' @param w vector with misclassified covariate
#'
#' @return a vector with probabilities
#'
#'
new_pi_conditional <- function(alpha, z, MC_0, MC_1, w){

  if(is.null(ncol(z)) || ncol(z) == 0){
    eta <- alpha[1]
  }else{
    eta <- alpha[1] + t(alpha[-c(1)]) %*% z
  }
  pi <- 1 / (1 + exp(-eta))

  #p(w=1) and p(w=0) that will be used as normalizing constants below
  M_22 <- MC_0[2,2]*(1-z) + MC_1[2,2]*z
  M_12 <- MC_0[1,2]*(1-z) + MC_1[1,2]*z
  pw1 <- M_22 * pi + M_12*(1 - pi)
  pw0 <- 1 - pw1

  # probabilities to sample x in the MC iterations
  sample_pi <- ifelse(w == 1, # if w=1
                      M_22*pi / pw1, # use this
                      (1 - M_22)*pi / pw0 # otherwise this
  )
  return(sample_pi)
}
