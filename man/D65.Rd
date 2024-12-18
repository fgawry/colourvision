\name{D65}
\alias{D65}
\docType{data}
\title{
CIE Standard Illuminant D65 in quantum flux (umol/m2/s)
}
\description{
CIE Standard Illuminant D65. Datum has already been converted to quantum flux (umol/m2/s) and therefore can be used in colour vision models directly.
}
\usage{data("D65")}
\format{
  A data frame with 107 observations on the following 2 variables.
  \describe{
    \item{\code{l.nm}}{a numeric vector}
    \item{\code{Standard.Illuminant.D65}}{a numeric vector}
  }
}
\source{
\url{https://www.cie.co.at}
}
\examples{
data("D65")
plot(D65, type="l")
}
\keyword{datasets}
