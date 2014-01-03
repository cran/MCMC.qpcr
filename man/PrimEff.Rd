\name{PrimEff}
\alias{PrimEff}
\title{
Determines qPCR amplification efficiencies from dilution series
}
\description{
Runs linear regression on Cq versus log2(RNA concentration), plots graph,
reports slope (ideally should be -1), and efficiency (with 95 percent confidence 
limits)
}
\usage{
PrimEff(data, plot = TRUE)
}
\arguments{
  \item{data}{
a dataframe containing three columns. First is RNA concentration. This could be
absolute as well as relative concentration (1/dilution factor). 
Second is the Cq value. Third is gene name. Replicate the same name across 
all the corresponding RNA concentrations. The dataframe may contain data for 
multiple genes.
}
  \item{plot}{
set plot=FALSE if the plot is not required
}
}
\details{
Run with at least 8 2-fold dilutions per gene
}
\value{
Plots the regression and under it, the values of slope and efficiency 
(plus and minus one SD). The dataframe may contain 
data for multiple genes, which will all be plotted 
together (so the reasonable limit is something like 25 genes)

Also returns a dataframe with columns: gene, efficiency, plus one SD, minus one SD, and intercept.
}
\references{
Matz MV, Wright RM, Scott JG (2013) No Control Genes Required: Bayesian Analysis of qRT-PCR Data. PLoS ONE 8(8): e71448. doi:10.1371/journal.pone.0071448
}
\author{
Mikhail V. Matz, UT Austin
<matz@utexas.edu>
}
\examples{

data(dilutions)
PrimEff(dilutions)

}