\name{getNormalizedData}
\alias{getNormalizedData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extracts qPCR model predictions
}
\description{
Generates a table of model-derived log2-transformed transcript abundances without global sample effects (i.e., corresponding to efficiency-corrected and normalized qPCR data)
}
\usage{
getNormalizedData(model, data, controls=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
qPCR model: the output of mcmc.qpcr or mcmc.qpcr.lognormal function fitted with two additional options: random="sample", pr=TRUE . These options do not change the inferences of main effects but make it possible to retain among-sample variation of expression for each gene while still subtracting the global sample effects (i.e., perform "normalization")
}
  \item{data}{
The dataset that was analysed to generate the model (output of cq2counts or cq2log functions)
}
  \item{controls}{
List of control genes; required if the mcmc.qpcr model was fit with the option normalize=TRUE
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
library(MCMC.qpcr)

# loading Cq data and amplification efficiencies
data(coral.stress) 
data(amp.eff) 

genecolumns=c(5,6,16,17) # specifying columns corresponding to genes of interest
conditions=c(1:4) # specifying columns containing factors  

# calculating molecule counts and reformatting:
dd=cq2counts(data=coral.stress,genecols=genecolumns,
condcols=conditions,effic=amp.eff,Cq1=37) 

# fitting the model (must include random="sample", pr=TRUE options)
mm=mcmc.qpcr(
	fixed="condition",
	data=dd,
	controls=c("nd5","rpl11"),
	nitt=4000,
	pr=TRUE,
	random="sample"
)

# extracting model predictions
pp=getNormalizedData(mm,dd)

# here is the normalized data:
pp$normData

# and here are the corresponding conditions:
pp$conditions

# putting them together for plotting:
ppcombo=cbind(stack(pp$normData),rep(pp$conditions))
names(ppcombo)[1:2]=c("expression","gene")

# plotting boxplots of normalized data:
ggplot(ppcombo,aes(condition,expression,colour=timepoint))+
	geom_boxplot()+
	facet_wrap(~gene,scales="free")+
	theme_bw()
}
