\name{MCMC.qpcr-package}
\alias{MCMC.qpcr}
\docType{package}
\title{
Bayesian analysis of qRT-PCR data
}
\description{
This package implements generalized linear mixed model analysis of qRT-PCR data 
so that the increase of variance towards higher Ct values is properly dealt with, 
and the lack of amplification is informative (function mcmc.qpcr). Sample-loading effects, 
gene-specific variances, and responses of all genes to each factor combination are all jointly 
estimated within a single model. The control genes can be specified as priors, with 
adjustable degree of expected stability. The analysis also works well without any 
control gene specifications. 

For higher-abundance datasets, a lognormal model is implemented that does not require 
Cq to counts conversion (function mcmc.qpcr.lognormal).

For higher-abundance datasets datasets in which the quality and/or quantity of RNA samples 
varies systematically (rather than randomly) across conditions, the analysis based on 
multigene normalization is implemented (function mcmc.qpcr.classic).

The package includes several functions for plotting the results and calculating statistical
significance (HPDplot, HPDplotBygene, HPDplotBygeneBygroup).

The detailed step-by-step tutorial is here:
http://www.bio.utexas.edu/research/matz_lab/matzlab/Methods_files/mcmc.qpcr.tutorial.release1.pdf
}
\details{
\tabular{ll}{
Package: mcmc.qpcr\cr
Type: Package\cr
Version: 1.0.1\cr
Date: 2013-07-23\cr
License: GPL-3\cr
}
}
\author{
Mikhail V. Matz, University of Texas at Austin
<matz@utexas.edu>
}
\references{
Matz MV, Wright RM, Scott JG (2013) No Control Genes Required: Bayesian Analysis of qRT-PCR Data. PLoS ONE 8(8): e71448. doi:10.1371/journal.pone.0071448
}
\keyword{ package }
\examples{

# here is the whole analysis using part of the coral stress dataset
# (use MCMC.qpcr tutorial for more extended examples)

# loading Cq data and amplification efficiencies
data(coral.stress) 
data(amp.eff) 
# extracting a subset of data 
cs.short=subset(coral.stress, timepoint=="one")

genecolumns=c(5,6,16,17) # specifying columns corresponding to genes of interest
conditions=c(1:4) # specifying columns containing factors  

# calculating molecule counts and reformatting:
dd=cq2counts(data=cs.short,genecols=genecolumns,
condcols=conditions,effic=amp.eff,Cq1=37) 

# fitting the model
mm=mcmc.qpcr(
	fixed="condition",
	data=dd,
	controls=c("nd5","rpl11"),
	nitt=4000 # remove this line when analyzing real data!
)

# plotting log2(fold change) in response to heat stress for all genes
HPDplot(model=mm,factors="conditionheat",main="response to heat stress")

# plotting abundances of individual genes across all conditions
# step 1: defining conditions
cds=list(
  control=list(factors=0), # gene-specific intercept
  stress=list(factors=c(0,"conditionheat")) # multiple effects will be summed up
  )

# step 2: plotting gene after gene on the same panel
HPDplotBygene(model=mm,gene="actin",conditions=cds,col="cyan3",
pch=17,jitter=-0.1,ylim=c(-3.5,15),pval="z")
HPDplotBygene(model=mm,gene="hsp16",conditions=cds,
newplot=FALSE,col="coral",pch=19,jitter=0.1,pval="z")

# step 3: adding legend
legend(0.5,10,"actin",lty=1,col="cyan3",pch=17,bty="n")
legend(0.5,7,"hsp16",lty=1,col="coral",pch=19,bty="n")

}
