\name{cq2log}
\alias{cq2log}
\title{
Prepares qRT-PCR data for mcmc.qpcr analysis using lognormal and "classic" (normalization-based) models
}
\description{
Log-transforms and efficiency-corrects Cq values, converting them to natural logs fo relative abundances,
and stacks the dataset
}
\usage{
cq2log(data, genecols, condcols, effic, noamp = 38)
}
\arguments{
  \item{data}{
Raw qRT-PCR dataset, one Cq column per gene, plus columns containing factors. The Cq 
columns, in addition to the proper Cq values, may contain NA (missing data) and -1, 
which means no amplification observed (i.e., zero target molecules at the start of qPCR 
reaction). 
Column headers are either gene names or factor names. 
Any number of fixed factors is allowed; any number of random factors that are gene-specific 
scalars (such as effect of genotype, or block)
Must have a column called "sample", denoting individual cDNA preps.
Technical replicates should not be averaged, they should be represented as independent
rows with the same sample ID.
}
  \item{genecols}{
columns containing Cq data
}
  \item{condcols}{
columns corresponding to factors, including "sample" factor
}
  \item{effic}{
The PCR efficiency data for each of the analyzed genes. This is data frame with two columns:
gene name (must exactly match the headers of gene columns in Cq data table!) and efficiency (fold-
amplification per PCR cycle, determined from qPCR of serial dilutions; see PrimEff() function )
}
  \item{noamp}{
Value to replace instances of no amplification with. These instances would be coded by -1 in the data
table. Specify 'noamp=NA' if you want to disregard them, but by default they will be converted into an arbitrarily low value, 38  
}
}
\details{
The models that process cq2log output are expected to work well for datasets that don't have too many Cq values above 30 and don't have instances of no amplification. For examples, see cq2counts function.	
}
\value{
Returns a dataframe with a single response variable column ("count", even though it is actually a log-transformed relative abundance value), gene column ("gene") and several columns containing factors.
}
\references{
Matz MV, Wright RM, Scott JG (2013) No Control Genes Required: Bayesian Analysis of qRT-PCR Data. PLoS ONE 8(8): e71448. doi:10.1371/journal.pone.0071448
}
\author{
Mikhail V. Matz, UT Austin
}
