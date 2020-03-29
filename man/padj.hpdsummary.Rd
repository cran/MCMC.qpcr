\name{padj.hpdsummary}
\alias{padj.hpdsummary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Adjusts p-values within an HPDsummary() object for multiple comparisons 
}
\description{
Replaces raw p-values in an object returned by HPDsummary() 
by adjusted p-values corrected for multiple comparisons.
Disregards the entries corresponding to control genes.
}
\usage{
padj.hpdsummary(hpdsumm, controls = NULL, method = "BH")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hpdsumm}{
Output of HPDsummary()
}
  \item{controls}{
A vector of control gene names
}
  \item{method}{
p-value correction method (see function p.adjust), default is Benjamini-Hochberg
}
}
\value{
HPDsummary object with original p-values replaced by the corrected ones.
}
\references{
Matz MV, Wright RM, Scott JG (2013) No Control Genes Required: 
Bayesian Analysis of qRT-PCR Data. PLoS ONE 8(8): e71448.
doi:10.1371/journal.pone.0071448
}
\author{
Mikhail V. Matz, matz@mail.utexas.edu
}
