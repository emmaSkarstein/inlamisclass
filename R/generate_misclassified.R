#' Simulate data
#'
#' @inheritParams generate_misclassified
#'
#' @return A data frame with the response, one continuous covariate and one discrete covariate, and the probabilities associated with the discrete covariate.
#' @keywords internal
generate_data <- function(n, p = 2, betas = c(1, 2, 2),
                               alphas = c(0.5, 0.2), family_y = "rnorm",
                               sd_y = 1){

  # For now just one perfectly measured covariate. Could be generalized to p-1.
  if(p-1 == 1){
    z <- stats::runif(n, -1, 1)
  }else{
    z <- rep(0, n)
  }

  # Generate binary covariate
  eta_x <- alphas[1] + alphas[2]*z
  pi_x <- 1/(1 + exp(-eta_x))
  x <- stats::rbinom(n = n, size = 1, prob = pi_x)

  # Generate response
  eta_y <- betas[1] + betas[2]*x + betas[3]*z
  if(family_y == "rnorm"){
    y <- stats::rnorm(n, mean = eta_y, sd = sd_y)
  }else{
    pi_y <- 1/(1 + exp(-eta_y))
    y <- stats::rbinom(n, 1, pi_y)
  }

  data <- data.frame(y = y, z = z, x = x, pi_x = pi_x)
  return(data)
}

#' Generate misclassification in one covariate
#'
#' @param data a data frame from the function `generate_data()`.
#' @inheritParams generate_misclassified
#'
#' @return A data frame that also contains a column with the misclassified binary covariate.
#' @keywords internal
misclassify <- function(data, MC_matrix){
  n <- nrow(data)

  # Generate misclassified variable w
  # Misclassification probability depends on value of x
  w <- rep(NA, n)
  for(i in 1:n){
    p_misclass <- ifelse(data$x[i] == 1, MC_matrix[2,2], MC_matrix[2,1])
    w[i] <- stats::rbinom(1, 1, p_misclass)
  }

  data$w <- w

  return(data)
}

#' Simulate data with a misclassified covariate
#'
#' @param n number of observations
#' @param p number of covariates
#' @param MC_matrix misclassification matrix
#' @param betas vector of coefficients for the model of interest (the order is 1. intercept, 2. coefficient for variable with error, 3. coefficient for variable without error)
#' @param alphas vector of coefficients for the imputation model (the order is 1. intercept, 2. coefficient for variable without error)
#' @param family_y likelihood family for the model of interest
#' @param sd_y standard deviation for the response variable
#'
#' @return A data frame with a response variable, a continuous covariate, a binary covariate with misclassification, the correct version of that binary covariate, and the probabilities associated with the binary covariate.
#' @export
#'
generate_misclassified <- function(n, p = 2,
                                   MC_matrix,
                                   betas = c(1, 2, 2),
                                   alphas = c(0.5, 0.2),
                                   family_y = "rnorm",
                                   sd_y = 1){
  data <- generate_data(n = n, p = p, betas = betas, alphas = alphas,
                        family_y = family_y, sd_y = sd_y)
  misclass_data <- misclassify(data, MC_matrix)
  return(misclass_data)
}


