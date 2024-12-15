\name{logistic}
\alias{logistic}
\title{
Logistic curve
}
\description{
Generates a logistic curve.}
\usage{
logistic(x = seq(300, 700, 1), x0, L, k)
}
\arguments{
  \item{x}{Range.
}
  \item{x0}{Midpoint value.
}
  \item{L}{Maximum value.
}
  \item{k}{Steepness of the curve
}
}
\references{
\url{https://en.wikipedia.org/wiki/Logistic_function}
}
\author{
Felipe M. Gawryszewski \email{f.gawry@gmail.com}
}
\examples{
l<-logistic(x=seq(300,700,1), x0=650, L=50, k=0.04)
plot(l, type="l")
}