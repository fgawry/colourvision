\name{bee}
\alias{bee}
\docType{data}
\title{
Honeybee photoreceptors
}
\description{
Honeybee (Apis mellifera) photoreceptor sensitivity curves.
}
\usage{data("bee")}
\format{
  A data frame with 401 observations on the following 4 variables.
  \describe{
    \item{\code{Wavelength}}{a numeric vector}
    \item{\code{UV}}{a numeric vector}
    \item{\code{Blue}}{a numeric vector}
    \item{\code{Green}}{a numeric vector}
  }
}
\details{
Original data were interpolated to 1nm intervals from 300 to 700nm.
}
\source{
Chittka, L., and P. Kevan. 2005. Flower colour as advertisement. Pp. 157-196 in Practical pollination biology.}
\examples{
data("bee")
plot(bee[,2]~bee[,1], col = "violet", type="l", xlab="Wavelength(nm)", ylab= "Absorbance")
lines(bee[,3]~bee[,1], col = "blue", type="l")
lines(bee[,4]~bee[,1], col = "green", type="l")
}
\keyword{datasets}