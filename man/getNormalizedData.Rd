\name{getNormalizedData}
\alias{getNormalizedData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extracts qPCR model predictions
}
\description{
Generates a table of model-derived log2-transformed transcript abundances without sample effects (i.e., corresponding to efficiency-corrected and normalized qPCR data)
}
\usage{
getNormalizedData(model, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
qPCR model: the output of mcmc.qpcr or mcmc.qpcr.lognormal function fitted with pr=TRUE option 
}
  \item{data}{
The dataset that was analysed to generate the model (output of cq2counts or cq2log functions)
}
}
\value{
The function returns a list of two data frames. The first one, normData, is the model-predicted log2-transformed transcript abundances table. It has one column per gene and one row per sample. The second data frame, conditions, is a table of experimental conditions corresponding to the normData table.
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

# extracting model predictions
pp=getNormalizedData(mm,dd)
# displaying model prediction
pp$normData
# displaying experimental conditions
pp$conditions

}
