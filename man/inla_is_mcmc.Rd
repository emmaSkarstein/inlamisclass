% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inla_mcmc.R
\name{inla_is_mcmc}
\alias{inla_is_mcmc}
\title{Fit model with misclassified covariate with INLA, importance sampling and MCMC}
\usage{
inla_is_mcmc(
  formula_moi,
  formula_imp = NULL,
  alpha0,
  MC_matrix,
  data,
  niter = 1600,
  ...
)
}
\arguments{
\item{formula_moi}{formula for the model of interest}

\item{formula_imp}{formula for the imputation model}

\item{alpha0}{initial values for coefficients for imputation model}

\item{MC_matrix}{misclassification matrix}

\item{data}{data for the model}

\item{niter}{number of iterations for the MCMC}

\item{...}{further arguments to be passed to inla().}
}
\value{
A list with the resulting model.
}
\description{
Fit model with misclassified covariate with INLA, importance sampling and MCMC
}
