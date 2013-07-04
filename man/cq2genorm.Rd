\name{cq2genorm}
\alias{cq2genorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reformats raw Ct data for geNorm analysis (non-parametric selection of
stable control genes) as implemented in selectHKgenes function (package SLqPCR)
}
\description{
Does similar procedures as cq2log, but in the end converts the data into
relative expression values, leaves it unstacked, and removes factor columns.
}
\usage{
cq2genorm(data, genes, effic, noamp=38)
}
\arguments{
  \item{data}{
raw qRT-PCR data; see help for cq2counts for details on formatting
}
  \item{genes}{
Vector of names of the potential control genes. Make sure you select only the  potential control genes for this analysis, otherwise geNorm might prefer actual responsive genes if they co-vary and have very low variance.
}
  \item{effic}{
PCR efficiency data for each of the analyzed genes; see help for cq2counts for details on formatting. 
}
  \item{noamp}{
what to do about Ct values '-1', denoting no amplification. By default, these will be assigned 
an arbitrarily low expression value corresponding to slightly less than one molecule. Specify 
'noamp=NA' if you want to skip samples containing any of these.
}
}
\details{
See cq2counts help page for details.
}
\value{
A dataset to be fed into the function selectHKgenes (package SLqPCR).
}
\references{
pending
}
\author{
Mikhail V. Matz, UT Austin
}
