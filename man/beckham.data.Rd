\name{beckham.data}
\alias{beckham.data}
\docType{data}
\title{
Cellular heat stress response data.
}
\description{
A typical qRT-PCR dataset where a series of treatments are compared to the global control.
}
\usage{data(beckham.data)}
\format{
  A data frame with 45 observations on the following 13 variables.
  \describe{
    \item{\code{sample}}{cDNA sample (biological replicate): a factor with levels \code{a} \code{b} \code{c} \code{d} \code{e} \code{f} \code{g} \code{h} \code{i} \code{j} \code{k} \code{l} \code{m} \code{n} \code{o}}
    \item{\code{tr}}{treatment: a factor with levels \code{ctl}:control \code{hs}:heat stress \code{pshs}:pre-stress followed by heat stress}
    \item{\code{time}}{hours of heat stress exposure: a numeric vector}
    \item{\code{gapdh}}{a numeric vector}
    \item{\code{hsp110}}{a numeric vector}
    \item{\code{hspb}}{a numeric vector}
    \item{\code{egr}}{a numeric vector}
    \item{\code{gadd}}{a numeric vector}
    \item{\code{dnajb1}}{a numeric vector}
    \item{\code{dnajb4}}{a numeric vector}
    \item{\code{atf}}{a numeric vector}
    \item{\code{dnaja4}}{a numeric vector}
    \item{\code{fos}}{a numeric vector}
  }
}
\details{
Includes global control (ctl), heat stressed cells (hs), and cells pre-stressed by mild heating prior to heat stress (pshs).
Two time points: one hour and three hours.
}
\source{
data for Figure 6 in Beckham et al. Microarray analysis of cellular thermotolerance. 
Lasers Surg Med. 2010 Dec;42(10):752-65. doi: 10.1002/lsm.20983
}
\examples{
data(beckham.data)
}
