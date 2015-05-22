\name{mcmc.qpcr.classic}
\alias{mcmc.qpcr.classic}
\title{
Analyzes qRT-PCR data using "classic" model, based on multigene normalization.
}
\description{
Normalizes the data using specified control genes, fits a single model to estimate changes at all genes. Use for datasets with not too many Cq values above 30.
}
\usage{
mcmc.qpcr.classic(fixed = NULL, globalFixed = NULL, random = NULL, globalRandom = NULL,
data, controls, genebysample = TRUE, geneSpecRes=FALSE, center = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fixed}{
desired combination of fixed effects, as a text string. Do not use "*" symbol, list
it fully, such as: 'factor1+factor2+factor1:factor2'. 
}
  \item{globalFixed}{
Vector of fixed covariates (categorical or continuous) that are expected to affect
all genes in the sample in the same way. These would be typically related to quality 
and/or quantity of RNA, such as RIN value.
}
  \item{random}{
A vector of names for gene-specific scalar random effects, such as 'c("effect1","effect2")'.
}
  \item{globalRandom}{
Random covariates (categorical only) affecting all genes, similar to globalFixed. 
}
  \item{data}{
output of the cq2log() function
}
  \item{controls}{
Vector of control gene names. These will be pushed to the back of the gene list during model fitting,
in the reverse order.
}
  \item{genebysample}{
Whether random gene by sample interactions should be
modeled as an additional random effect. If the model fails to converge, specify 'genebysample=F'. 
}
  \item{geneSpecRes}{
Whether the model should include gene-specific residuals.
This was the default in MCMC.qpcr v.1.0; now it is switched off
since it does not have any visible effect on the results (as long as genebysample=TRUE), 
and only makes the model converge slower or fail to converge.  
}
  \item{center}{
Whether to center the normalized log-transformed Cq values arond the mean for each gene. Centering does not affect the inference; it only makes the plots of the results more comprehensible (in my opinion). 
}
 \item{\dots}{
other options for MCMCglmm function, such as nitt (number of iterations), thin (tinning interval), and burnin (number of initial iterations to disregard). For a more precise inference (but longer runs) specify 'nitt=45000, thin=20, burnin=5000'. See MCMCglmm documentation for more details.
}
}
\details{
This function takes an as input the log-transformed relative abundance values, performs multigene normalization as per Vandesompele et al 2002, and runs an MCMC chain to fit a lognormal linear mixed model to estimate gene expression changes jointly at all genes. It is very powerful, as long as the average stability of control genes can be trusted. 

This function requires data prepared by cq2log function and must have the control genes specified; otherwise the arguments and syntax are similar to the mcmc.qpcr function.
}
\value{
An MCMCglmm object. See mcmc.qpcr function for details and examples.
}
\references{
Matz MV, Wright RM, Scott JG (2013) No Control Genes Required: Bayesian Analysis of qRT-PCR Data. PLoS ONE 8(8): e71448. doi:10.1371/journal.pone.0071448
}
\author{
Mikhail V. Matz, University of Texas at Austin
}

