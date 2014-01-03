\name{HPDpoints}
\alias{HPDpoints}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
HPDplot, HPDpoints
}
\description{
Calculates and plots posterior means with 95\% confidence intervals for specified 
fixed effects (or their combination) for all genes. HPDpoints only adds graphs to an existing plot.
}
\usage{
HPDpoints(model, factors, factors2 = NULL, ylimits = NULL, 
hpdtype = "w", inverse = F, jitter = 0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
The output of mcmc.qpcr function.
}
  \item{factors}{
A vector of names of fixed effects of interest; see details.
}
  \item{factors2}{
A second vector of fixed effect names to be subtracted from the first; see details.
}
  \item{ylimits}{
Y-limits for the plot such as c(-3,6); autoscale by default.
}
  \item{hpdtype}{
Specify hpdtype="l" to plot the upper and lower 95\% confidence limits as a continuous
dashed line across all genes. By default (hpdtype="w") the limits are plotted as whiskers around
each point.
}
  \item{inverse}{
Plot the inverse of the result.
}
  \item{jitter}{
For hpdtype="w", shifts the plotted values and whiskers by the specified distance along
the x axis (reasonable jitter values are 0.15 or -0.15, for example). 
}
  \item{\dots}{
Various plot() options; such as col (color of lines and symbols), pch (type of symbol),
main (plot title) etc. 
}
}
\details{
See details in HPDplot()
}
\value{
A graph added to a plot.
}
\references{
Matz MV, Wright RM, Scott JG (2013) No Control Genes Required: Bayesian Analysis of qRT-PCR Data. PLoS ONE 8(8): e71448. doi:10.1371/journal.pone.0071448
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
	nitt=4000 # remove this line when analyzing real data!
)

# plotting log2(fold change) in response to heat stress for all genes
HPDplot(model=mm,factors="conditionheat",main="response to heat stress")

}