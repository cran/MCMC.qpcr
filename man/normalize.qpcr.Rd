\name{normalize.qpcr}
\alias{normalize.qpcr}
\title{
Internal function called by mcmc.qpcr.classic
}
\description{
Performs multigene normalization as per Vandesompele et al 2002, centers the resulting gene expression values around the gene's means
}
\usage{
normalize.qpcr(data, controls, center = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
The output of cq2log()
}
  \item{controls}{
Vector of control genes, such as 'controls=c("eif3h","nd5","rpl11")'
}
  \item{center}{
Whether to perform centering around mean or not.
}
}
\value{
normalized dataset
}
\references{
pending
}
\author{
Mikhail Matz, UT Austin
}

