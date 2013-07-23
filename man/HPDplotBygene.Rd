\name{HPDplotBygene}
\alias{HPDplotBygene}
\title{
Plots qPCR analysis results for individual genes.
}
\description{
For a particular gene, plots model-predicted values (and confidence intervals) for a series of specified  fixed effect combinations ("conditions").
}
\usage{
HPDplotBygene(model, gene, conditions, pval = "mcmc", newplot = T, 
ylimits = NULL, inverse = F, jitter = 0, plot = T, yscale = "log2", 
interval = "ci", grid = F, zero = F, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
model object produced by MCMC.qpcr()
}
  \item{gene}{
name of the gene to plot
}
  \item{conditions}{
A list naming the conditions to plot and defining them as combination of fixed effects (See example below). '0' denotes gene-specific intercept.
}
  \item{pval}{
Type of p-value to report: 'mcmc' - MCMC-based (default), 'z' - based on Bayesian z-score. Use 'z' to approximate p-values lower than 2/[MCMC sample size] (with default MCMC.qpcr settings, this comes to 0.002) 
}
  \item{newplot}{
When TRUE, a new plot should be created. When FALSE, or a graph will be added to an existing plot.
}
  \item{ylimits}{
Y-limits for the plot such as c(-3,6); autoscale by default.
}
  \item{inverse}{
When TRUE, the inverse of the data will be plotted.
}
  \item{jitter}{
Shifts the plotted values and whiskers by the specified distance along
the x axis (reasonable jitter values are 0.15 or -0.15, for example). This helps plot 
several graphs on the same plot without overlapping. 
}
  \item{plot}{
When FALSE, no plot will be generated; the function will just list mean pairwise differences and p-values.
}
  \item{yscale}{
Scale on which to represent the data. In all mcmc.qpcr models the model scale is natural logarithm, which
I prefer to translate into log2 or log10 (if the differences are orders of magnitude) for better human readability. 
The default is 'log2'; other options are 'log10' and 'native' (no rescaling of the model data). 
There is also a beta-option 'proportion', which is not useful for qPCR. It was added to cannibalize HPDplotBygene 
function for plotting results of the model operating with arcsin-square root
transfromed proportions. With yscale="proportions", the plot will be on the original proportion scale but the tukey-like 
differences will still be reported on the asin(sqrt()) transformed scale.
}
  \item{interval}{
'ci' (default) will plot 95\% confidence limits of the posterior distribution, 'sd' will plot the mean plus/minus one standard deviation of the posterior.
}
  \item{grid}{
When TRUE, a vertical grid separating conditions will be added
}
  \item{zero}{
When TRUE, a y=0 line will be added.
}
  \item{\dots}{
Various plot() options.
}
}
\value{
Generates a point-whiskers plot, lists pairwise mean differenes between all conditions, calculates and lists pairwise p-values (not corrected for multiple testing). 
}
\references{
pending
}
\author{
Mikhal V. Matz, UT Austin, matz@utexas.edu
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