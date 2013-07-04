\name{HPDplotBygeneBygroup}
\alias{HPDplotBygeneBygroup}
\title{
Plots qPCR analysis results for individual genes
}
\description{
For a specified gene, makes overlayed plots such as produced by HPDplotBygene()
}
\usage{
HPDplotBygeneBygroup(model, gene, group1, group2, group3 = NULL, 
interval = "ci", colors = c("coral", "cyan3", "grey50"), 
symbols = c(19, 17, 15), jitter = 0.16, yscale = "log2", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
model object produced by mcmc.qpcr()
}
  \item{gene}{
name of the gene to plot
}
  \item{group1}{
Combination of factors defining the first group (see HPDplotBygene() for details). 
}
  \item{group2}{
Combination of factors defining the second group.
}
  \item{group3}{
(optional) Combination of factors defining the third group.
}
  \item{interval}{
'ci' (default) will plot 95\% confidence limits of the posterior distribution, 'sd' will plot the mean plus/minus one standard deviation of the posterior.
}
  \item{colors}{
Colors to use for different groups (see ?par -> col).
}
  \item{symbols}{
Symbols to use for different groups (see ?par -> pch).
}
  \item{jitter}{
Jitter distance between groups.
}
  \item{yscale}{
For qPCR analysis, the default 'log2' should be used. Other options are 'log10', 'proportion' (assumes that arcsin-suqare root transformation was used and plots the results on the original scale), and 'absolute' (no data transformation for plotting).  
}
  \item{\dots}{
additional parameters for HPDplotBygene() function, such as pval (see HPDplotBygene() help)
}
}
\value{
Generates a point-whiskers plot, lists pairwise mean differenes between all conditions, calculates and lists pairwise p-values (not corrected for multiple testing). 
}
\references{
pending
}
\author{
Mikhal V. Matz, UT Austin
<matz@utexas.edu>
}

