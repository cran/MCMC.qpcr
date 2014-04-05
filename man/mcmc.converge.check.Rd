\name{mcmc.converge.check}
\alias{mcmc.converge.check}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MCMC diagnostic plots
}
\description{
A wrapper function for plot.MCMCglmm to plot diagnostic convergence plots for selected fixed effects
}
\usage{
mcmc.converge.check(model, factors, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
output of mcmc.qpcr (or any MCMCglmm class object)
}
  \item{factors}{
A vector of names of fixed effects of interest; see details in HPDplot help page.
}
  \item{\dots}{
other options to pass to plot.MCMCglmm
}
}
\value{
A series of plots for each gene-specific fixed effect.

The MCMC trace plot is on the left, to see if there is convergence (lack of systematic
trend) and no autocorrelation (no low-frequency waves). If lack of convergence is suspected, 
try increasing number of iterations and burnin by specifying, for example, nitt=50000,
burnin=5000, as additional options for mcmc.qpcr. If autocorrelation is present, increase 
thinning interval by specifying thin=20 in mcmc.qpcr (you might wish to increase the number of iterations,
nitt, to keep the size of MCMC sample the same)

The right plot is posterior density distribution.
}
\references{
Matz MV, Wright RM, Scott JG (2013) No Control Genes Required: Bayesian Analysis of qRT-PCR Data. PLoS ONE 8(8): e71448. doi:10.1371/journal.pone.0071448
}
\author{
Mikhail V. Matz, UT Austin
}
