\name{outlierSamples}
\alias{outlierSamples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
detects outlier samples in qPCR data 
}
\description{
reports samples that have too little starting material relative to others (by default, less by two standard deviations) 
}
\usage{
outlierSamples(model, data, z.cutoff = -2)
}
\arguments{
  \item{model}{
qPCR model: the output of mcmc.qpcr or mcmc.qpcr.lognormal function fitted with pr=TRUE option 
}
  \item{data}{
The dataset that was analysed to generate the model (output of cq2counts or cq2log functions)
}
  \item{z.cutoff}{
z-score cutoff to report an outlier sample.
}
}
\value{
A vector containing outlier sample names. 
}
\references{
Matz MV, Wright RM, Scott JG (2013) No Control Genes Required: Bayesian Analysis of qRT-PCR Data. PLoS ONE 8(8): e71448. doi:10.1371/journal.pone.0071448
}
\author{
Mikhail V. Matz, University of Texas at Austin
<matz@utexas.edu>
}
\examples{
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
	nitt=4000, # remove this line when analyzing real data!
	pr=TRUE
)

# detecting outliers
outliers=outlierSamples(mm,dd)

# removing outliers
dd=dd[!(dd$sample \%in\% outliers),]

}
