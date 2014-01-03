\name{padj.qpcr}
\alias{padj.qpcr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates adjusted p-values corrected for multiple comparisons
}
\description{
Takes the output of HPDplot(...,plot=FALSE), adds two columns of adjusted p-values (z and mcmc)
Disregards the entries corresponding to control genes
}
\usage{
padj.qpcr(data, controls = NULL, method = "BH")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Output of HPDplot(...,plot=FALSE); may be several several such outputs concatenated with rbind
}
  \item{controls}{
A vector of control gene names
}
  \item{method}{
p-value correction method (see function p.adjust), default is Benjamini-Hochberg
}
}
\value{
The dataframe derived from the input, with added columns "padj.z" and "padj.mcmc".
See tutorial for examples.
}
\references{
Matz MV, Wright RM, Scott JG (2013) No Control Genes Required: Bayesian Analysis of qRT-PCR Data. PLoS ONE 8(8): e71448. doi:10.1371/journal.pone.0071448
}
\author{
Mikhail V. Matz, matz@mail.utexas.edu
}
