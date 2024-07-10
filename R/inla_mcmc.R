#' Fit model with misclassified covariate with INLA and MCMC
#'
#' @param formula_moi formula for the model of interest
#' @param formula_imp formula for the imputation model
#' @param alpha initial values for coefficients for imputation model
#' @param MC_matrix misclassification matrix
#' @param data data for the model
#' @param niter number of iterations for the MCMC
#' @param nburnin number of burn-in iterations
#'
#' @return A data frame with the resulting model.
#' @export
#'
inla_mcmc <- function(formula_moi, formula_imp = NULL,
                      alpha, MC_matrix,
                      data, niter = 1600, nburnin = 100){
  # So far, this assumes the exposure model just has one covariate.
  # And the formula_moi does not really do much either

  r.out.naive <- INLA::inla(formula_moi, data = data)

  models_list <- list()
  models_list[[1]] <- r.out.naive

  nn <- nrow(data)

  mc.out <- list()

  # Identify response variable
  response <- all.vars(formula_moi)[1]

  # Identify error variable
  error_var <- all.vars(formula_imp)[1]

  # Identify error free covariates
  error_free_covs <- labels(terms(formula_imp))


  for (ii in 1:niter){
    #print(sprintf("%d", ii))

    sample_pi <- new_pi(alpha, z = data[, error_free_covs],
                        MC_matrix, w = data[, error_var])

    # We have a model for x; use the sample_pi probabilities derived above to sample from it
    xstar <- stats::rbinom(nn, 1, sample_pi)

    new_data <- cbind(data, xstar = xstar)
    new_formula_moi <- stats::reformulate(response = response, termlabels = c("xstar", error_free_covs))

    #r.inla <- inla(y ~ xstar + z, data=dd,num.threads=1, control.mode = list(result = r.out,restart=TRUE))
    if(ii < niter){
      r.inla <- INLA::inla(new_formula_moi,
                     data = new_data,
                     control.mode = list(result = models_list[[ii]],
                                         restart = TRUE),
                     control.compute = list(return.marginals = FALSE),
                     control.inla = list(int.strategy = 'eb'))
    }else if(ii == niter){
      r.inla <- INLA::inla(new_formula_moi,
                     data = new_data,
                     control.mode = list(result = models_list[[ii]],
                                         restart = TRUE))
    }

    models_list[[ii+1]] <- r.inla

    # sample new alpha for the exposure model of x. To this end, use xstar to fit the model and then extract alpha
    new_formula_imp <- stats::reformulate(response = "xstar", termlabels = error_free_covs)
    r.inla.x <- INLA::inla(new_formula_imp, data = new_data,
                     family = "binomial",
                     control.compute = list(config = TRUE,
                                            return.marginals = FALSE),
                     control.inla = list(int.strategy = 'eb'))

    r.alpha <- INLA::inla.posterior.sample(n = 1, r.inla.x)
    alpha <- r.alpha[[1]]$latent[c(nn+1, nn+2)]

    # results <- c(r.inla$mlik[1,1] + sum(log(lik_error)), r.inla$summary.fixed$mean)
    results <- c(r.inla$mlik[1,1],
                 r.inla.x$mlik[1,1],
                 r.inla$summary.fixed$mean,
                 r.inla$summary.fixed$`0.025quant`[2],
                 r.inla$summary.fixed$`0.975quant`[2],
                 alpha[1],
                 alpha[2])
    names(results) <- c("mlik", "mlik_z", "intercept", "betax", "betaz",
                        "betax0.025", "betax0.975", "alpha0", "alphaz")

    mc.out[[ii]] <- results
  }

  #mc.results <- list(mcmc_results = mc.out,
  #                   naive_model = summary(lm(y ~ w + z, data = data))$coef,
  #                   correct_model = summary(lm(y ~ x + z, data = data))$coef,
  #                   exposure_model = summary(glm(x ~ z, family = "binomial", data = data)))

  #save(mc.results,file="code/RData/mc.out.uncertainty.RData")
  return(mc.out)
}

#' Calculation of new probabilities for x, given the alpha (and rest)
#'
#' @param alpha coefficients for the imputation model
#' @param z vector (if only one covariate) or matrix of covariates for imputation model.
#' @param MC_matrix misclassification matrix
#' @param w vector with misclassified covariate
#'
#' @return a vector with probabilities
#' @export
#'
new_pi <- function(alpha, z, MC_matrix, w){

  eta <- alpha[1] + t(alpha[-c(1)]) %*% z
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
