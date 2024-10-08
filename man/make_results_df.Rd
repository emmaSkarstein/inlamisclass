% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary_and_plot.R
\name{make_results_df}
\alias{make_results_df}
\title{Summarize results from inla_is or inla_mcmc}
\usage{
make_results_df(
  mcmc_results,
  niter = length(mcmc_results),
  nburnin = length(mcmc_results) - niter
)
}
\arguments{
\item{mcmc_results}{results from inla_is or inla_mcmc}

\item{niter}{number of iterations for the MCMC}

\item{nburnin}{number of iterations for the burn-in}
}
\value{
A dataframe with summary statistics.
}
\description{
Summarize results from inla_is or inla_mcmc
}
