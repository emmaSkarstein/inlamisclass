#' Fit model with misclassified covariate with INLA and importance sampling
#'
#' @param formula_moi formula for the model of interest
#' @param formula_imp formula for the imputation model
#' @param MC_matrix misclassification matrix
#' @param data data for the model
#' @param niter number of iterations for the MCMC
#' @param alpha values for coefficients of imputation model
#' @param ncores number of cores for the parallel computation.
#' @param conditional if the misclassification probabilities are conditional on some variable, then this should be the name of that variable.
#' @param missing_only TRUE/FALSE indicating that there is no misclassification, only missingness in a covariate. By default FALSE. Note that if there is misclassification AND missingness, this should still be FALSE.
#' @param ... further arguments to be passed to inla().
#'
#' @return A data frame with the resulting model.
#' @export
#'
inla_is_misclass <- function(formula_moi, formula_imp = NULL,
                    alpha, MC_matrix = NULL, data, niter, ncores = 4,
                    conditional = NULL, missing_only = FALSE, ...){
  if(is.null(MC_matrix) && !missing_only){
    stop("You have not provided a misclassification matrix. Do you mean to do missing data imputation only? If so, specify the argument 'missing_only = TRUE'.")
  }
  if(!is.null(MC_matrix) && missing_only){
    stop("You have provided a misclassification matrix even though you have also specified 'missing_only = TRUE', meaning that there should not be any misclassification.")
  }
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
    # If missing_only == TRUE, we want to use the true covariate values whenever the covariate is not missing.
    if(missing_only){
      missing_ind <- is.na(data[, vars$error_var])
      xstar[!missing_ind] <- data[, vars$error_var]
    }

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

#' Identify roles of covariates in formula
#'
#' @inheritParams inla_is_misclass
#'
#' @return A list of named variables.
#'
#' @keywords internal
identify_vars <- function(formula_moi, formula_imp){
  # Identify response variable
  response <- all.vars(formula_moi)[1]

  # Identify error variable
  error_var <- all.vars(formula_imp)[1]

  # Covariates in imputation model
  imp_covs <- labels(stats::terms(formula_imp))

  # Identify error free covariates
  error_free_covs <- setdiff(all.vars(formula_moi), c(response, error_var))

  return(list(response = response,
              error_var = error_var,
              imp_covs = imp_covs,
              error_free_covs = error_free_covs))
}

#' Calculate success probability for binary variable depending on covariate(s)
#'
#' @inheritParams new_pi
#'
#' @return a probability vector
#' @keywords internal

calculate_pi <- function(alpha, z){
  # Need to check nrow(t(z)) since
  # - if z is a vector, ncol(z) = NULL, but nrow(t(z)) = 1
  # - if z is a matrix, nrow(t(z)) = ncol(z)
  # - if z is empty (no covariates in imp. model), then nrow(t(z)) = 0.
  if(nrow(t(z)) == 0){
    eta <- alpha[1]
  }else{
    eta <- alpha[1] + t(alpha[-c(1)]) %*% z
  }
  pi <- 1 / (1 + exp(-eta))

  return(pi)
}

#' Calculation of new probabilities for x, given the alpha (and rest)
#'
#' @param alpha coefficients for the imputation model
#' @param z vector (if only one covariate) or matrix of covariates for imputation model.
#' @param MC_matrix misclassification matrix
#' @param w vector with misclassified covariate
#'
#' @return a vector with probabilities
#' @keywords internal
#'
new_pi <- function(alpha, z, MC_matrix, w){
  # pi = success probability of x
  pi <- calculate_pi(alpha, z)

  # If MC_matrix is NULL, assume this is missing data imputation. Then, just return pi.
  if(is.null(MC_matrix)){
    return(pi)
  }

  pw1 <- MC_matrix[1, 2]*(1 - pi) + MC_matrix[2, 2] * pi
  pw0 <- MC_matrix[1, 1]*(1 - pi) + MC_matrix[2, 1] * pi

  # Probabilities to sample x in the MC iterations
  sample_pi <- ifelse(w == 1, # if w=1
                      MC_matrix[2, 2]*pi / pw1, # use this
                      (1 - MC_matrix[2, 2])*pi / pw0 # otherwise this
  )

  #return(data.frame(sample_pi = sample_pi, sample_pi_old = sample_pi_old, pw1_old = pw1_old[1,], pw1 = pw1[1,], pw0_old = pw0_old[1,], pw0 = pw0[1,]))
  return(sample_pi = sample_pi)
}


#' Calculation of new probabilities for x, given the alpha (and rest)
#'
#' @param alpha coefficients for the imputation model
#' @param z vector of covariates for imputation model, and that the misclassification matrix is conditioned on.
#' @param MC_matrix list of two misclassification matrices, one for conditioning_var = 0, and the other for conditioning_var = 1.
#' @param w vector with misclassified covariate
#' @param conditioning_var variable that the misclassification is conditional on
#'
#' @return a vector with probabilities
#'
#'
new_pi_conditional <- function(alpha, z, MC_matrix, w, conditioning_var){
  # pi = success probability of x
  pi <- calculate_pi(alpha, z)

  # If MC_matrix is NULL, assume this is missing data imputation. Then, just return pi.
  if(is.null(MC_matrix)){
    return(pi)
  }

  # Extract the two MC matrices
  MC_0 <- MC_matrix$MC_0
  MC_1 <- MC_matrix$MC_1

  # This produces a list with n entries,
  M_22 <- MC_0[2,2]*(1-conditioning_var) + MC_1[2,2]*conditioning_var
  M_12 <- MC_0[2,1]*(1-conditioning_var) + MC_1[2,1]*conditioning_var

  #p(w=1) and p(w=0) that will be used as normalizing constants below
  pw1 <- M_22 * pi + M_12*(1 - pi)
  pw0 <- 1 - pw1

  # probabilities to sample x in the MC iterations
  sample_pi <- ifelse(w == 1, # if w=1
                      M_22*pi / pw1, # use this
                      (1 - M_22)*pi / pw0 # otherwise this
  )
  return(sample_pi)
}

