\name{EMline}
\alias{EMline}
\title{
Endler and Mielke (2005) 1-D colour space
}
\description{
Plots a colour space for dichromatic Endler and Mielke (2005) colour vision model.
}
\usage{
EMline(x,y=rep(0, length(x)), type="length",
         vnames=c("E1","E2"), 
         ylim="auto", xlim="auto",
         ann=FALSE, axes = FALSE, ...)
}
\arguments{
  \item{x}{x-coordinate}
  \item{y}{y-coordinate}
  \item{type}{Whether the colour space should be built with a fixed vector length, or a fixed edge length (distance between vertices).}
  \item{vnames}{Vector names}
  \item{ylim}{see \code{\link{plot}} function.}
  \item{xlim}{see \code{\link{plot}} function.}
  \item{ann}{see \code{\link{par}} function.}
  \item{axes}{see \code{\link{plot.default}} function.}
  \item{...}{Other arguments passed to \code{\link{plot}}.}   
}
\details{
The original model is available for tetrachromatic animals only. Colour space is built either with a vector length = 0.75 or a edge length = sqrt(3/2), to match the tetrahedron proposed by Endler and Mielke (2005).
}
\references{
Endler, J. A., and P. Mielke. 2005. Comparing entire colour patterns as birds see them. Biol J Linn Soc 86:405-431.
}

\author{
Felipe M. Gawryszewski \email{f.gawry@gmail.com}
}
\seealso{
\code{\link{EMtriangle}}, \code{\link{EMtetrahedron}},  \code{\link{EMmodel}}
}
\examples{
EMline(x=0.1, pch=16, col="red")
}