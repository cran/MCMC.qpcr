\name{softNorm}
\alias{softNorm}
\title{
Accessory function to mcmc.qpcr() to perform soft normalization
}
\description{
adds a fake NORM gene with counts equal to geometric mean of conrtol genes,
removes control genes from the dataset
}
\usage{
softNorm(data,controls)
}
\arguments{
  \item{data}{
The dataset; output of cq2counts()
}
  \item{controls}{
Vector fo control gene names
}
}
\value{
A modified dataset with NORM gene set as the reference, and no control genes
}
\author{
Mikhail Matz, UT Austin
<matz@utexas.edu>
}
