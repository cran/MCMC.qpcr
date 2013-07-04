\name{dilutions}
\alias{dilutions}
\docType{data}
\title{
Data to determine amplification efficiency
}
\description{
Cq data for a series of four-fold dilutions for two targets; the input for PrimEff()
}
\usage{data(dilutions)}
\format{
  A data frame with 76 observations on the following 3 variables.
  \describe{
    \item{\code{dna}}{a numeric vector}
    \item{\code{cq}}{a numeric vector}
    \item{\code{gene}}{a factor with levels \code{chrom} \code{eif3h}}
  }
}
\examples{
data(dilutions)
PrimEff(dilutions)
}
\keyword{datasets}
