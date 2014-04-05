\name{HPDplot}
\alias{HPDplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting fixed effects for all genes for a single combination of factors
}
\description{
Calculates and plots posterior means with 95\% credible intervals for specified 
fixed effects (or their combination) for all genes
}
\usage{
HPDplot(model, factors, factors2 = NULL, ylimits = NULL, 
hpdtype = "w", inverse = F, jitter = 0, plot = T, grid = T, zero = T, ...)
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
Specify hpdtype="l" to plot the upper and lower 95\% credible limits as a continuous
dashed line across all genes. This is useful to compare credible intervals among several models on the same plot. 
By default (hpdtype="w") the limits are plotted as whiskers around each point.
}
  \item{inverse}{
Plot the inverse of the result.
}
  \item{jitter}{
For hpdtype="w", shifts the plotted values and whiskers by the specified distance along
the x axis (reasonable jitter values are 0.15 or -0.15, for example). This helps plot 
several results (different models or factor combinations) on the same plot (use 
HPDpoints to add to existing plot) 
}
  \item{plot}{
if plot = FALSE the function returns a table of calculated posterior modes, means, upper 
and lower 95\% credible limits (all on log(2) scale), and two types of p-values: derived from Bayesian z-scores, 
and derived directly from MCMC sample. 

All such outputs for a given experiment should be concatenated
with rbind and processed by padj.qpcr() function to adjust the p-values for multiple comparisons 
(disregarding the entries corresponding to control genes)
}
  \item{grid}{
Whether to draw vertical grid lines to separate genes. 
}
  \item{zero}{
Whether to draw a horizontal line at 0. 
}
  \item{\dots}{
Various plot() options; such as col (color of lines and symbols), pch (type of symbol),
main (plot title) etc. 
}
}
\details{
Use summary(MCMCglmm object) first to see what fixed effect names are actually used 
in the output. For example, if summary shows:

gene1:conditionheat\cr
gene2:conditionheat\cr
....\cr
gene1:timepointtwo\cr
gene2:timepointtwo\cr
....\cr
gene1:conditionheat:timepointtwo\cr
gene2:conditionheat:timepointtwo

, it is possible to specify factors="conditionheat" to plot only the effects of the heat.

If a vector of several fixed effect names is given, for example: 
factors=c("timepointtwo","treatmentheat:timepointtwo")
the function will plot the posterior mean and credible interval for the sum of these 
effects.

If a second vector is also given, for example,\cr 
factors=c("f1","f2"), factors2=c("f3","f4")\cr
the function will plot the difference between the sums of these two groups of factors.
This is useful for pairwise analysis of differences in complicated designs.
}
\value{
A plot or a table (plot = F).

Use the function HPDpoints() if you need to add graphs to already existing plot.

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
	nitt=3000,burnin=2000 # remove this line when analyzing real data!
)

# plotting log2(fold change) in response to heat stress for all genes
HPDplot(model=mm,factors="conditionheat",main="response to heat stress")

}
