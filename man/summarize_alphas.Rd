% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary_and_plot.R
\name{summarize_alphas}
\alias{summarize_alphas}
\title{Calculate summary statistics for the alphas}
\usage{
summarize_alphas(mcmc_results, niter = NULL, nburnin)
}
\arguments{
\item{mcmc_results}{results from inla_is or inla_mcmc}

\item{niter}{number of iterations for the MCMC}

\item{nburnin}{number of iterations for the burn-in}
}
\value{
Summary statistics for the alphas
}
\description{
Calculate summary statistics for the alphas
}
