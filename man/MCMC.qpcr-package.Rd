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
http://www.bio.utexas.edu/research/matz_lab/matzlab/Methods_files/mcmc.qpcr.tutorial.pdf
}
\details{
\tabular{ll}{
Package: MCMC.qpcr\cr
Type: Package\cr
Version: 1.2\cr
Date: 2015-05-22\cr
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
