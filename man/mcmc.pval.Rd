\name{mcmc.pval}
\alias{mcmc.pval}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
calculates p-value based on Bayesian z-score
}
\description{
Accessory function to HPDplot and varHPDplot series
}
\usage{
mcmc.pval(dat, testlim = 0, sided = 2)
}
\arguments{
  \item{dat}{
a table of MCMC samples
}
  \item{testlim}{
positive value to test whether the posterior crosses it (for variance components)
}
  \item{sided}{
sidedness of the test
}
}
\value{
A vector of calculated Bayesian p-values. The minimal possible p-value is 1/(MCMC sample size)
}
\author{
Mikhail Matz, UT Austin
<matz@utexas.edu>
}
