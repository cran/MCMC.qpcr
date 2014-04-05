\name{mcmc.qpcr}
\alias{mcmc.qpcr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Analyzes qRT-PCR data using generalized linear mixed model
}
\description{
Wrapper function for MCMCglmm by Jarrod Hadfield, designed for qRT-PCR data 
analysis. 
}
\usage{
mcmc.qpcr(fixed, random = NULL, data, controls = NULL, normalize=FALSE,include = NULL, 
m.fix = 1.2, v.fix = NULL, geneSpecRes=TRUE, vprior="uninf", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fixed}{
desired combination of fixed effects, as a text string. Do not use "*" symbol, list
it fully, such as: 'factor1+factor2+factor1:factor2'. 
}
  \item{random}{
A vector of names for gene-specific scalar random effects, such as 'c("effect1","effect2")'.
}
  \item{data}{
output of the cq2counts() function
}
  \item{controls}{
Vector of control gene names. These will be pushed to the back of the gene list during model fitting,
in the reverse order. Controls are NOT NECESSARY for the analysis.
}
  \item{normalize}{
If controls are specified, requiring 'normalize=TRUE' will perform "soft normalization":
the geometric mean of control genes will be used to infer global effects 
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
other options for MCMCglmm function, such as nitt (number of iterations), thin (tinning interval), and burnin (number of initial iterations to disregard). For a more precise inference specify 'nitt=45000, thin=20, burnin=5000'. See MCMCglmm documentation for more details.
}
}
\details{
This function constructs priors and runs an MCMC chain to fit a Poisson-lognormal generalized linear mixed model to the counts data.

The fixed effects for the model by default assume a gene-specific intercept, and gene-specific effect for each of the listed fixed factors. 

The user-specified random effects are all assumed to be gene-specific with no covariances. 

There is also a universal random factor: the scalar random effect of sample,
which accounts for the unequal template loading. 

Residual variances are assumed to be gene-specific with no covariances, with weakly informative inverse Wishart prior (variance=1, nu=[number of genes]-0.998).

The priors for fixed effects are diffuse gaussians with 
a mean at 0 and very large variances (1e+8), except for control genes, for which the prior variances are defined 
by the m.fix parameter. For the gene-specific random effects and residual variation, non-informative priors are used to achieve results equivalent to ML estimation. For control genes, when v.fix parameter is specified, it will be used to restrict residial variance. 

Both m.fix and v.fix parameters should be specified as allowed average fold-change, so the lowest they can go is 1 (no variation). 

All control genes share the same m.fix and v.fix parameters.
}
\value{
An MCMCglmm object. HPDsummary() function within this package summaizes these data, 
calculates all gene-wise credible intervals and p-values, and plots the results either 
as line-point-whiskers graph or a bar-whiskers graph using ggplot2 functions. 

HPDsummary() only works for experiments with a single multilevel factor or two fully 
crossed multilevel factors. Use finctions HPDplot(), HPDplotBygene() and HPDplotBygeneBygroup() 
to summarize and plot more complicated designs.

For more useful operations on MCMCglmm objects, such as posterior.mode(), HPDinterval(),
and plot(), see documentation for MCMCglmm package.
}
\references{
Matz MV, Wright RM, Scott JG (2013) No Control Genes Required: Bayesian Analysis of qRT-PCR Data. PLoS ONE 8(8): e71448. doi:10.1371/journal.pone.0071448
}
\author{
Mikhail V. Matz, University of Texas at Austin
<matz@utexas.edu>
}
\examples{

data(beckham.data)
data(beckham.eff)

# analysing the first 5 genes 
# (to try it with all 10 genes, change the line below to gcol=4:13)
gcol=4:8 
ccol=1:3 # columns containing experimental conditions

# recalculating into molecule counts, reformatting
qs=cq2counts(data=beckham.data,genecols=gcol,
condcols=ccol,effic=beckham.eff,Cq1=37)

# creating a single factor, 'treatment.time', out of 'tr' and 'time'
qs$treatment.time=as.factor(paste(qs$tr,qs$time,sep="."))

# fitting a naive model
naive=mcmc.qpcr(
	fixed="treatment.time",
	data=qs,
	nitt=3000,burnin=2000 # remove this line in actual analysis!
)

#summary plot of inferred abundances
#s1=HPDsummary(model=naive,data=qs)

#summary plot of fold-changes relative to the global control
s0=HPDsummary(model=naive,data=qs,relative=TRUE)

# pairwise differences and their significances for each gene:
s0$geneWise

}


