\name{coral.stress}
\alias{coral.stress}
\docType{data}
\title{
RT-qPCR of stress response in coral Porites astreoides
}
\description{
Timepoint one is one day of heat-light stress; timepoint two is recovery from it
the next day. Controls are the fragments never exposed to stress. Individual denotes
a coral colony that yielded several fragments for different treatments. 
}
\usage{data(coral.stress)}
\format{
  A data frame with 64 observations on the following 19 variables.
  \describe{
    \item{\code{sample}}{a numeric vector}
    \item{\code{individual}}{a factor with levels \code{s1} \code{s12} \code{s13} \code{s15} \code{s2} \code{s4} \code{s5} \code{s8}}
    \item{\code{condition}}{a factor with levels \code{control} \code{heat}}
    \item{\code{timepoint}}{a factor with levels \code{one} \code{two}}
    \item{\code{hsp16}}{a numeric vector}
    \item{\code{actin}}{a numeric vector}
    \item{\code{adk}}{a numeric vector}
    \item{\code{c3}}{a numeric vector}
    \item{\code{chrom}}{a numeric vector}
    \item{\code{clect}}{a numeric vector}
    \item{\code{eif3h}}{a numeric vector}
    \item{\code{g3pdh}}{a numeric vector}
    \item{\code{gsp2}}{a numeric vector}
    \item{\code{hsp60}}{a numeric vector}
    \item{\code{hsp90}}{a numeric vector}
    \item{\code{nd5}}{a numeric vector}
    \item{\code{rpl11}}{a numeric vector}
    \item{\code{spon2}}{a numeric vector}
    \item{\code{ubl3}}{a numeric vector}
  }
}
\examples{
data(coral.stress)
}
\keyword{datasets}
