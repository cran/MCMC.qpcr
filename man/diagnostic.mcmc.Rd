\name{diagnostic.mcmc}
\alias{diagnostic.mcmc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots three diagnostic plots to check the validity of the assumptions of linear model analysis.
}
\description{
Predicted vs observed plot tests for linearity, Scale-location plot tests for homoscedasticity, and Normal QQ plot tests for normality of the residuals.
} 
\usage{
diagnostic.mcmc(model, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
MCMCglmm object (a model fitted by mcmc.qpcr or mcmc.qpcr.gauss), obtained with additional options, 'pl=T, pr=T'
}
  \item{\dots}{
Various plot() options to modify color, shape and size of the plotteed points.
}
}
\value{
A plot with three panels.	
}
\references{
pending
}
\author{
Mikhail V. Matz, UT Austin
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
	pr=TRUE,pl=TRUE, # these flags are necessary for diagnostics
	nitt=4000 # remove this line when analyzing real data!
)
diagnostic.mcmc(mm)

}
