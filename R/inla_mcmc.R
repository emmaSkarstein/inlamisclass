#' Fpi
#'
#' @inheritParams Fm
#' @return pi (what is pi?)
#'
#' @keywords internal
Fpi <- function(X, beta, n){
  x <- cbind(rep(1, n), X)
  eta <- x %*% beta
  pi <- exp(eta)/(1+exp(eta))
  return(pi)
}

#' FyTilde
#'
#' @inheritParams Fm
#' @param Ntrials Number of trials (for when?)
#' @return yTilde (what is yTilde?)
#'
#' @keywords internal
FyTilde <- function(X, beta, Y, Ntrials = 1, n){
  pi <- Fpi(X = X, beta = beta, n = n)
  x <- cbind(rep(1, n), X)
  eta <- x %*% beta
  yTilde <- eta + (Y/Ntrials-pi)*1/(pi*(1-pi))
  return(yTilde)
}

#' FWMat
#'
#' @inheritParams Fm
#' @return WMat (what is WMat?)
#'
#' @keywords internal
FWMat <- function(X, beta, n){
  pi <- Fpi(X = X, beta = beta, n = n)
  WiInv <- 1/(pi*(1-pi))
  WMat <- diag(n)
  diag(WMat) <- 1/WiInv
  return(WMat) # Faktor zum Tunen der Acceptance Rate (Ueberlegen, warum man das darf)
}

#' FC
#'
#' @inheritParams Fm
#' @return C (what is C?)
#'
#' @keywords internal
FC <- function(X, beta, R, a, varB, n, p){
  WMat <- FWMat(X = X, beta = beta, n = n)
  x <- cbind(rep(1, n), X)
  C <- solve(1/varB*diag(1, p+1) + t(x) %*% WMat %*% x)
  return(C)
}

#' Fm
#'
#' @param X Covariate matrix
#' @param beta Coefficient vector
#' @param R not sure
#' @param a not sure
#' @param Y Proposed response value
#' @param varB Variance of what?
#'
#' @return m (what is m?)
#'
#' @keywords internal
Fm <- function(X, beta, R, a, Y, varB, n, p){
  WMat <- FWMat(X, beta, n)
  yTilde <- FyTilde(X = X, beta = beta, Y = Y, n = n)
  C <- FC(X = X, beta = beta, R = R, a = a, varB = varB, n = n, p = p)
  x <- cbind(rep(1, n), X)
  m <- t(C %*% (1/varB*a + t(x) %*% WMat %*% yTilde))
  return(m)
}

#' Fit model with misclassified covariate with INLA, importance sampling and MCMC
#'
#' @param formula_moi formula for the model of interest
#' @param formula_imp formula for the imputation model
#' @param alpha0 initial values for coefficients for imputation model
#' @param MC_matrix misclassification matrix
#' @param data data for the model
#' @param niter number of iterations for the MCMC
#' @param ... further arguments to be passed to inla().
#'
#' @return A list with the resulting model.
#' @export
#'
inla_is_mcmc <- function(formula_moi, formula_imp = NULL,
                         alpha0, MC_matrix,
                         data, niter = 1600, ...){

  r.out.naive <- INLA::inla(formula_moi, data = data)#, ...)

  alpha <- alpha0

  vars <- identify_vars(formula_moi = formula_moi, formula_imp = formula_imp)
  z <- data[, vars$imp_covs]
  w <- data[, vars$error_var]

  # Number of variables in the regression model for x
  p <- length(vars$imp_covs)

  n <- nrow(data)

  # Prior for alpha
  a <- rep(0, p+1)
  R <- diag(p+1)
  varAlpha <- 10
  diag(R) <- rep(varAlpha, p+1)

  Bb <- 1 #0.25 # Tune the acceptance rate


  mc.out <- list()
  for (ii in 1:niter){
    sample_pi <- new_pi(alpha, z = z,
                        MC_matrix, w = w)
    # We have a model for x; use the sample_pi probabilities derived above to sample from it
    xstar <- stats::rbinom(n, 1, sample_pi)

    new_data <- cbind(data, xstar = xstar)
    new_formula_moi <- stats::reformulate(response = vars$response, termlabels = c("xstar", vars$error_free_covs))

    #r.inla <- inla(y ~ xstar + z, data=dd,num.threads=1, control.mode = list(result = r.out,restart=TRUE))
    r.inla <- INLA::inla(new_formula_moi,
                   data=new_data,
                   control.mode = list(result = r.out.naive, restart=TRUE),
                   control.compute = list(return.marginals = FALSE),
                   control.inla = list(int.strategy = 'eb'))

    # Sample new alpha for the exposure model of x:
    C <- FC(X = z, beta = alpha, R = R, a = a, varB = varAlpha, n = n, p = p)
    m <- Fm(X = z, beta = alpha, R = R, a = a, Y = xstar, varB = varAlpha, n = n, p = p)
    pi <- Fpi(X = z, beta = alpha, n = n)

    alphaStar <- MASS::mvrnorm(n = 1, mu = m, Sigma = C*Bb)

    CStar <- FC(X = z, beta = alphaStar, R = R, a = a, varB = varAlpha, n = n, p = p)
    mStar <- Fm(X = z, beta = alphaStar, R = R, a = a, Y = xstar, varB = varAlpha, n = n, p = p)
    piStar <- Fpi(X = z, beta = alphaStar, n = n)

    likelyStar <- sum(stats::dbinom(xstar, 1, piStar, log=T))
    likely <- sum(stats::dbinom(xstar, 1, pi, log=T))

    logAlphaStar <- -1/2 * (alphaStar - a) %*% solve(R) %*% (alphaStar - a) + likelyStar

    logAlpha <- -1/2 * (alpha - a) %*% solve(R) %*% (alpha - a) + likely

    logQAlpha <- mvtnorm::dmvnorm(x = alpha, mean = mStar, sigma = CStar*Bb, log=T)

    logQAlphaStar <- mvtnorm::dmvnorm(x = alphaStar, mean = m, sigma = C*Bb, log=T)

    logacc <- logAlphaStar - logAlpha + logQAlpha - logQAlphaStar
    alphayes <- 0

    if(log(stats::runif(1)) <= logacc){
      alpha <- alphaStar
      alphayes <- 1
    }

    summary.fixed <- r.inla$summary.fixed
    variables <- rownames(summary.fixed)

    # Change name of xstar to error variable
    error_var_index <- which(variables == "xstar")
    variables[error_var_index] <- vars$error_var

    summary.fixed$variable <- variables

    # Set names to alpha vector
    alpha_names <- c("alpha.0", paste0("alpha.", vars$imp_covs))
    names(alpha) <- alpha_names

    results <- list(mlik = r.inla$mlik[1,1],
                    summary.fixed = summary.fixed,
                    alpha = alpha,
                    alphayes = alphayes,
                    accprob = exp(logacc))

    mc.out[[ii]] <- results
  }

  return(mc.out)
}
