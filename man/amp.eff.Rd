\name{amp.eff}
\alias{amp.eff}
\docType{data}
\title{
amplification efficiencies and experimental Cq1 (optional column)
}
\description{
Fold-amplification within each PCR cycle for each qRT-PCR assay. Determined by qPCR of
serial dilutions of the template (see Pfaffl 2001, Nucleic Acids Res 29:00)
}
\usage{data(amp.eff)}
\format{
  A data frame with 17 observations on the following 3 variables.
  \describe{
    \item{\code{gene}}{a factor with levels \code{actin} \code{adk} \code{c3} \code{chrom} \code{clect} \code{eif3h} \code{g3pdh} \code{gsp2} \code{hsp16} \code{hsp60} \code{hsp90} \code{nd5} \code{r18s} \code{rpl11} \code{spon2} \code{tgoln} \code{ubl3}}
    \item{\code{efficiency}}{a numeric vector}
    \item{\code{Cq1}}{a numeric vector}
  }
}
\examples{
data(amp.eff)
}
\keyword{datasets}
