#' Identify roles of covariates in formula
#'
#' @inheritParams inla_is_mcmc
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
  #pw1_old <- MC_matrix[2, 2] * pi + MC_matrix[1, 2]*(1 - pi)
  #pw0_old <- 1 - pw1_old

  #sample_pi_old <- ifelse(w == 1, # if w=1
  #                    MC_matrix[2, 2]*pi / pw1_old, # use this
  #                    (1 - MC_matrix[2, 2])*pi / pw0_old # otherwise this
  #)

  # Tentative:
  #pw1 <- MC_matrix[1, 1]*(1 - pi) + MC_matrix[1, 2] * pi
  #pw0 <- MC_matrix[2, 1]*(1 - pi) + MC_matrix[2, 2] * pi

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

#' Calculate some charateristics from simulated data
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
