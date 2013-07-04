\name{varHPDplot}
\alias{varHPDplot}
\title{
Plotting variance components for all genes... Use with caution: variance components are very sensitive to priors!
}
\description{
Calculates and plots posterior means with 95\% confidence intervals for specified 
variance components (or their combination, such as heritability) for all genes
}
\usage{
varHPDplot(model, factors = "units", factors2 = NULL, ylimits = NULL, 
hpdtype = "w", jitter = 0, plot = T, testlimit = 0.05, grid = T, zero = T, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
The output of mcmc.qpcr function.
}
  \item{factors}{
A vector of names of the variance components to be summed and plotted; default is residual variance.
}
  \item{factors2}{
A second vector of variance component names. If it is specified, the function will plot the proportion of 
variance attributable to the sum of these components relative to total variance ('factors2' plus 'factors')
}
  \item{ylimits}{
Y-limits for the plot, such as c(0,1); autoscale by default.
}
  \item{hpdtype}{
Specify hpdtype="l" to plot the upper and lower 95\% confidence limits as a continuous
dashed line across all genes. This is useful to evaluate the overall pattern of gene 
expression change. By default (hpdtype="w") the limits are plotted as whiskers around
each point.
}
  \item{jitter}{
For hpdtype="w", shifts the plotted values and whiskers by the specified distance along
the x axis (reasonable jitter values are 0.2 or -0.2, for example). This helps plot 
several results (different models or factor combinations) on the same plot (use 
HPDpoints to add to existing plot) 
}
  \item{plot}{
if plot = FALSE the function returns a table of calculated posterior means, upper 
and lower 95\% confidence limits of the deviations (square root of variances) on log(2) scale,   
and MCMC-based p-values to test whether the parameter exceeds some positive value specified as 'testlimit'. 
Proportions (when 'factors2' are specified) are based directly on variances.
}
  \item{testlimit}{
The lower boundary for testing whether the deviations (or proportions of variance) is significant
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
The variance components can be gene-specific residual variances ("units") or gene-specific random effects specified in the model.
Variance exploration might be useful for deciding which of the analyzed genes are the most stable. Also, proportion of variance attributable
to random effects associated with ancestry or individual identity can be used to infer heritabilities of gene expression. 
}
\value{
A plot or a table (plot = F). 

Use the function points.var.hpd() if you need to add graphs to already existing plot.
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

# fitting the model without control genes, with random effect of individual
mm0=mcmc.qpcr(
	fixed="condition",
	random="individual",
	data=dd,
#	thin=40, # un-remark this line for real data analysis
#	burnin=10000, # un-remark this line for real data analysis
	nitt=4000 # set this to 100000 when analyzing real data
)
# plotting broad-sense heritabilities (proportion of variance due to individual)
varHPDplot(factors="units",factors2="individual",model=mm0)

# does broad-sense heritability exceed some not-too-small value, say 0.1 ?
vv=varHPDplot(factors="units",factors2="individual",model=mm0,
plot=FALSE, testlimit=0.1)
vv

}