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
mcmc.qpcr(fixed, random = NULL, data, controls = NULL, include = NULL, 
m.fix = 1.2, v.fix = NULL, genebysample = T, vprior="flat", ...)
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
output of the cq2counts() function
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
Allowed residual fold-change deviation of the control genes between samples. If the dataset includes techincal replicates, 
this will apply to the gene by sample interaction term, if not, to the residual variation.
}
  \item{genebysample}{
Whether random gene by sample interactions should be
modeled as an additional random effect (see details). This requires technical replicates; if there are too many samples for which 
technical replicates are missing the model may not run, in which case specify 'genebysample=F'. 
}
   \item{vprior}{
Allows to set weakly informative inverse Wishart priors for variance components instead of the default flat prior. vprior="iw" will assume large prior variance (1), vprior="iw01" will assume small prior variance (0.1). Be careful with these, setting informative priors may affect the results if the data is not extensive enough.
}
  \item{\dots}{
other options for MCMCglmm function, such as nitt (number of iterations), thin (tinning interval), and burnin (number of initial iterations to disregard). For a more precise inference, such as required for model comparison, specify 'nitt=45000, thin=20, burnin=5000' (such a run may take
10-15 minutes, though). See MCMCglmm documentation for more details.
}
}
\details{
This function constructs priors and runs an MCMC chain to fit a Poisson-lognormal generalized linear mixed model to the counts data.

The fixed effects for the model by default assume a gene-specific intercept, and gene-specific effect for each of the listed fixed factors. 

The specified random effects are all assumed to be gene-specific with no covariances. 

There are also two universal random factors. The first one the scalar random effect of sample,
which accounts for the unequal template loading. The second one is random gene by sample interaction, 
to estimate how much variation is indeed due to difference between samples rather than to measurement error (the deviations between technical replicates). Note that the dataset must contain technical replicates for this variance component to work; 
if there are no replicates the factor will not be included irrespective of what the function call says. If only a few samples have technical
replicates, the model may not run; specify 'genebysample=F'.

Residual variances are assumed to be gene-specific.

The prior is set up to reflect the above structure. The priors for fixed effects are diffuse gaussians with 
a mean at 0 and very large variances (1e+8), except for control genes, for which the prior variances are defined 
by the m.fix parameter. For the gene-specific random effects and residual variation, non-informative priors are used to achieve results equivalent to ML estimation. For control genes, when v.fix parameter is specified, it will be used to restrict the variance due to gene by sample interaction
if the dataset contains technical replicates, and if it does not, it will restrict the residual variance. 

Both m.fix and v.fix parameters should be specified as allowed average fold-change, so the lowest they can go is 1 (no variation). 
All control genes share the same m.fix and v.fix parameters.
}
\value{
An MCMCglmm object. It contains all the results
of MCMC sampling during model fitting. If the object is named, for example, mm, use summary(mm) to see means, 
95\% credible intervals, and Bayesian p-values of all fixed effects, and how variance partitions 
between random effects and residuals. 
Note that the effects and credible intervals in the 
summary will be expressed as a natural logarithm of fold change; divide them by log(2) 
to convert to logarithm with the base 2 (at least for me, qPCR results are easier to understand on the log2 scale, since such values are like PCR cycles). 

HPDplot() function within this package plots the results for all genes across a particular condition and outputs p-values. HPDplotBygene() function plots results for a number of specified conditions for a particular gene and performs pairwise comparisons between conditions. Function HPDplotBygeneBygroup() overlays several graphs produced by HPDplotBygene() on the same plot.

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
	nitt=4000 # remove this line in actual analysis!
)

#summary plot of inferred abundances
#s1=HPDsummary(model=naive,data=qs)

#summary plot of fold-changes relative to the global control
s0=HPDsummary(model=naive,data=qs,relative=TRUE)

# pairwise differences and their significances for each gene:
s0$geneWise

}


