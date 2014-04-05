\name{mcmc.qpcr.lognormal}
\alias{mcmc.qpcr.lognormal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fits a lognormal linear mixed model to qRT-PCR data. 
}
\description{
Use in combination with cq2log(), on datasets without too many Cq values above 30. 
}
\usage{
mcmc.qpcr.lognormal(fixed, random = NULL, data, controls = NULL, 
normalize=FALSE, include = NULL, m.fix = 1.2, v.fix = NULL, 
geneSpecRes=TRUE,vprior = "uninf", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fixed}{
desired combination of fixed effects, as a text string. Do not use "*" symbol, list
it fully, such as: 'factor1+factor2+factor1:factor2'. 
}
  \item{random}{
A vector of names for gene-specific random effects, such as 'c("effect1","effect2")'.
}
  \item{data}{
output of the cq2log() function
}
  \item{controls}{
Vector of control gene names. These will be pushed to the back of the gene list during model fitting,
in the reverse order.
}
  \item{normalize}{
If controls are specified, requiring 'normalize=TRUE' will perform "soft normalization":
the mean of control genes will be used to infer global effects 
(common to all genes) across conditions; these will be subtracted when summarizing the data
with HPDsummary(). Use this option when template abundances are correlated with 
conditions. Note that this is different from normalizing data within each sample, as per
Vandesompele et al 2002; this would be "hard normalization" (use mcmc.qpcr.classic() if you want it).
}
  \item{include}{
How many of the control genes ('controls') to actually incorporate as priors during model
fitting. If left unspecified, all the 'controls' will be used. If 'include=0', the model
will be fitted without using any of the control genes as priors. If 'include' equals some number 
less than the number of 'controls', 
only the first 'include' of them will be used as priors. In all these cases, 
the 'controls' will appear in the same order in the output, in the end of the gene list
rather than according to their alphabetical position among other genes. This is useful when visually comparing the results
of models fitted with different number of control genes, using HPDplot and HPDpoints functions. 
}
  \item{m.fix}{
Allowed average fold-change of the control genes in response to any fixed factor combination.
}
  \item{v.fix}{
Allowed residual fold-change deviation of the control genes. Applies to the residual variation term.
}
  \item{geneSpecRes}{
Whether the model should include gene-specific residuals.
Keep it TRUE unless the model fails to converge.
}
   \item{vprior}{
The default prior is an uninformative inverse Wishart with assumed variance (V) at 1 
and the degree of belief parameter (nu) at 0. With 'prior="iw"' and 'prior="iw01"' nu is equal  [number of genes]-0.998, resulting in a weakly informative prior that is commonly used in this type of inference. vprior="iw" will assume large prior variance (1), vprior="iw01" will assume small prior variance (0.1). 
}
 \item{\dots}{
other options for MCMCglmm function, such as nitt (number of iterations), thin (tinning interval), and burnin (number of initial iterations to disregard). For a more precise inference (but longer runs) specify 'nitt=45000, thin=20, burnin=5000'. See MCMCglmm documentation for more details.
}
}
\details{
This function constructs priors and runs an MCMC chain to fit a lognormal linear mixed model to the log-transformed relative abundances data. The data for this function must be prepared by cq2log instead of cq2counts function; otherwise, the arguments and syntax are identical to mcmc.qpcr function.
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


