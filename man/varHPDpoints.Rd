\name{varHPDpoints}
\alias{varHPDpoints}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting variance components for all genes... Use with caution: variance components are very sensitive to priors!
}
\description{
Same as varHPDplot but adds graphs to an existing plot
}
\usage{
varHPDpoints(model, factors = "units", factors2 = NULL, hpdtype = "w", jitter = 0, ...)
}
\arguments{
  \item{model}{
The output of mcmc.qpcr function.
}
  \item{factors}{
A vector of names of the variance components to be summed and plotted; default is residual variance.
}
  \item{factors2}{
A second vector of variance component names. If it is specified, the function will plot the proportion of 
variance attributable to the sum of these components relative to total variance ('factors2' plus 'factors')
}
  \item{hpdtype}{
Specify hpdtype="l" to plot the upper and lower 95\% confidence limits as a continuous
dashed line across all genes. This is useful to evaluate the overall pattern of gene 
expression change. By default (hpdtype="w") the limits are plotted as whiskers around
each point.
}
  \item{jitter}{
For hpdtype="w", shifts the plotted values and whiskers by the specified distance along
the x axis (reasonable jitter values are 0.2 or -0.2, for example). This helps plot 
several results (different models or factor combinations) on the same plot (use 
HPDpoints to add to existing plot) 
}
  \item{\dots}{
Various plot() options; such as col (color of lines and symbols), pch (type of symbol),
main (plot title) etc. 
}
}
\details{
The variance components can be gene-specific residual variances ("units") or gene-specific random effects specified in the model.
Variance exploration might be useful for deciding which of the analyzed genes are the most stable. Also, proportion of variance attributable
to random effects associated with ancestry or individual identity can be used to infer heritabilities of gene expression.
See varHPDplot() for examples. 
}
\value{
A graph added to an existing plot. Use varHPDplot to generate a new plot, or produce a table of results (with plot=FALSE).
}
\references{
Matz MV, Wright RM, Scott JG (2013) No Control Genes Required: Bayesian Analysis of qRT-PCR Data. PLoS ONE 8(8): e71448. doi:10.1371/journal.pone.0071448
}
\author{
Mikhail V. Matz, UT Austin
<matz@utexas.edu>
}
