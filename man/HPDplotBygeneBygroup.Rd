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
Scale on which to represent the data. In all mcmc.qpcr models the model scale is natural logarithm, which
I prefer to translate into log2 or log10 (if the differences are orders of magnitude) for better human readability. 
The default is 'log2'; other options are 'log10' and 'native' (no rescaling of the model data). 
There is also a beta-option 'proportion', which is not useful for qPCR. It was added to cannibalize HPDplotBygene 
function for plotting results of the model operating with arcsin-square root
transfromed proportions. With yscale="proportions", the plot will be on the original proportion scale but the tukey-like 
differences will still be reported on the asin(sqrt()) transformed scale.
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

