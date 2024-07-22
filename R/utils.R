#' Identify roles of covariates in formula
#'
#' @inheritParams inla_mcmc
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
              errror_free_covs = error_free_covs))
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

  #p(w=1) and p(w=0) that will be used as normalizing constants below
  pw1 <- MC_matrix[2, 2] * pi + MC_matrix[1, 2]*(1 - pi)
  pw0 <- 1 - pw1

  # probabilities to sample x in the MC iterations
  sample_pi <- ifelse(w == 1, # if w=1
                      MC_matrix[2, 2]*pi / pw1, # use this
                      (1 - MC_matrix[2, 2])*pi / pw0 # otherwise this
  )
  return(sample_pi)
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
new_pi_conditional <- function(alpha, z, MC_matrix, w){

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

  MC_0 <- MC_matrix$MC_0
  MC_1 <- MC_matrix$MC_1

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
