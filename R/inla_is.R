#' Fit model with misclassified covariate with INLA and importance sampling
#'
#' @param data data for model
#' @param niter number of iterations to run
#' @param sample_pi probabilities to sample the misclassified covariate from.
#'
#' @return A data frame with the resulting model.
#' @export
#'
inla_is <- function(data, niter, sample_pi){
  r.out <- INLA::inla(y ~ x + z, data = data)
  r.out.naive <- INLA::inla(y ~ w + z, data = data)

  models_list <- list()
  models_list[[1]] <- r.out.naive

  n <- nrow(data)

  mc.out <- parallel::mclapply(1:niter, function(i) {

    # We have a model for x; use the sample_pi probabilities derived above to sample from it
    xstar <- stats::rbinom(n, 1, sample_pi)

    dd <- list(y = data$y, xstar = xstar, z = data$z)
    if(i < niter){
      r.inla <- INLA::inla(y ~ xstar + z,
                     data = dd,
                     num.threads = 1,
                     control.mode = list(result = r.out.naive,
                                         restart = TRUE),
                     control.compute = list(return.marginals = FALSE),
                     control.inla = list(int.strategy = 'eb')
      )
    }else if(i == niter){
      r.inla <- INLA::inla(y ~ xstar + z,
                     data = dd,
                     num.threads = 1,
                     control.mode = list(result = r.out.naive,
                                         restart = TRUE))
    }

    results <- c(r.inla$mlik[1,1],
                 r.inla$summary.fixed$mean,
                 r.inla$summary.fixed$`0.025quant`[2],
                 r.inla$summary.fixed$`0.975quant`[2])

    names(results) <- c("mlik","intercept","betax","betaz","betax0.025","betax0.975")

    results

  }, mc.cores = 4)
}
