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
include = NULL, m.fix = 1.2, v.fix = NULL, genebysample = T, vprior = "flat", ...)
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
Allowed residual fold-change deviation of the control genes between samples.
}
  \item{genebysample}{
Whether random gene by sample interactions should be
modeled as an additional random effect (see details). This requires technical replicates; if there are too many samples for which 
technical replicates are missing the model may not run, in which case specify 'genebysample=F'. 
}
   \item{vprior}{
Allows to set weakly informative inverse Wishart priors for variance components instead of the default flat prior. vprior="iw" will assume large prior variance (1), vprior="iw01" will assume small prior variance (0.1). Be careful with these, setting informative priors may affect the results if the data is not extensive enough, or all the variances are very small.
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
pending
}
\author{
Mikhail V. Matz, University of Texas at Austin
}


