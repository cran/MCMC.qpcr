\name{trellisByGene}
\alias{trellisByGene}
\title{
For two-way designs, plots mcmc.qpcr model predictions gene by gene
}
\description{
For each gene, plots model-predicted values and 95\% credible intervals.
}
\usage{
trellisByGene(modelSummary,xFactor,groupFactor,
nrow=2,lineWidth=0.4,whiskerWidth=0.2,pointSize=2.5,
facetScales="free_y",ylab="log(abundance)",
legendPos="bottom",posDodge=0.3)
}
\arguments{
  \item{modelSummary}{
two-way design model summary produced by HPDsummary()
}
  \item{xFactor}{
factor to form the x-axis
}
  \item{groupFactor}{
factor to form separate lines on the plot
}
  \item{nrow}{
number of rows in the resulting trellis plot
}
  \item{lineWidth}{
line width, passed as 'lwd' to geom_errorbar function (ggplot2)
}
  \item{whiskerWidth}{
width of the line denoting 95\% CI margin, passed as 'width' to geom_errorbar function (ggplot2)
}
  \item{pointSize}{
passed as 'size' to geom_point function of ggplot2
}
  \item{facetScales}{
passed as 'scales' to facet_wrap function of ggplot2
}
  \item{ylab}{
y-axis label
}
  \item{legendPos}{
passed as 'legend.position' to theme function of ggplot2
}
  \item{posDodge}{
position dodge, increase for more jitter
}
}
\value{
A ggplot2 type object
}
\references{
Matz MV, Wright RM, Scott JG (2013) No Control Genes Required: Bayesian Analysis of qRT-PCR Data. PLoS ONE 8(8): e71448. doi:10.1371/journal.pone.0071448
}
\author{
Mikhal V. Matz, UT Austin, matz@utexas.edu
}
\examples{

# loading Cq data and amplification efficiencies
data(coral.stress) 
data(amp.eff) 

genecolumns=c(5,6,16,17) # specifying columns corresponding to genes of interest
conditions=c(1:4) # specifying columns containing factors  

# calculating molecule counts and reformatting:
dd=cq2counts(data=coral.stress,genecols=genecolumns,
condcols=conditions,effic=amp.eff,Cq1=37) 

# fitting the 2-way model
mm=mcmc.qpcr(
	fixed="condition+timepoint+condition:timepoint",
	data=dd,
	nitt=4000 # remark this line to analyze real data!
)

# summarizing results
ss=HPDsummary(mm,data=dd,summ.plot=FALSE)

# plotting predicted means and 95\% CIs gene by gene
trellisByGene(ss,xFactor="condition",groupFactor="timepoint")
}