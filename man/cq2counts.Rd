\name{cq2counts}
\alias{cq2counts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prepares qRT-PCR data for mcmc.qpcr analysis
}
\description{
Converts Cq values into molecule counts,
and stacks the dataset
}
\usage{
cq2counts(data, genecols, condcols, effic, Cq1 = NULL)
}
%- maybe also 'usage' for other objects documented here.
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
columns that contain Cq data
}
  \item{condcols}{
columns corresponding to factors, including "sample" factor
}
  \item{effic}{
The PCR efficiency data for each of the analyzed genes. This is data frame with two columns:
gene name (must exactly match the headers of gene columns in Cq data table!) and efficiency (fold-
amplification per PCR cycle, determined from qPCR of serial dilutions; see PrimEff() function )
}
  \item{Cq1}{
The Cq of a single molecule. If left unspecified, it will be calculated from the efficiency (E)
using approximate formula Cq1=51.6-7.56*E, derived empirically for Roche's LightCycler 480. 
Cq1 does not seem to have much effect on relative quantification results unless it is 
wildly off (by 2-3 cycles). For an unknown qPCR instrument a single Cq1=37 could be assumed
for all genes. 

Note: If all experimental Cq values are less than 30, Cq1 variation (within a reasonable range, 35-39)
will not have any effect on the results whatsoever, so just go for Cq1=37.  
}
}
\value{
Returns a dataframe with a single response variable column ("count"), gene column 
("gene") and several columns containing factors.

Note: The purpose of Cq to counts conversion is to enable generalized linear model analysis,
which would take care of the heteroscedasticity and occasional 'empty' amplification trials 
for low-abundant targets. Although this works well, the absolute values of molecule counts returned by 
cq2counts are still approximate, so they should not be relied upon for true single-molecule analysis.
}
\references{
pending
}
\author{
Mikhail V. Matz, UT Austin
<matz@utexas.edu>
}
\examples{

data(coral.stress) 
data(amp.eff) 
genecolumns=c(5:19) # specifying where the Ct data are in the data table
conditions=c(1:4) # specifying data table columns containing factors  
# calculating molecule counts and reformatting:
dd=cq2counts(data=coral.stress,genecols=genecolumns,
condcols=conditions,effic=amp.eff,Cq1=37) 

}
