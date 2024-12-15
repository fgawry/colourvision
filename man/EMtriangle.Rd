\name{EMtriangle}
\alias{EMtriangle}
\title{
Endler and Mielke (2005) triangle colour space
}
\description{
Plots a triangle colour space for trichromatic Endler and Mielke (2005) colour vision model.
}
\usage{
EMtriangle(x, y, type=c("length", "edge"), vnames=c("u","s","m"),
           ylim=c(-0.9,0.9), xlim=c(-0.9,0.9),
           pch=16, bty="n",yaxt="n",xaxt="n",
           col="black", asp=1, ann=FALSE, vectors=FALSE, ...)
}
\arguments{
  \item{x}{x-coordinate}
  \item{y}{y-coordinate}
  \item{type}{Whether the colour space should be built with a fixed vector length, or a fixed edge length (distance between vertices).}
  \item{vnames}{Vector names}
  \item{ylim}{see \code{\link{plot}} function.}
  \item{xlim}{see \code{\link{plot}} function.}
  \item{pch}{see \code{\link{par}} function.}
  \item{bty}{see \code{\link{par}} function.}
  \item{yaxt}{see \code{\link{par}} function.}
  \item{xaxt}{see \code{\link{par}} function.}
  \item{col}{see \code{\link{par}} function.}
  \item{asp}{see \code{\link{plot}} function.}
  \item{ann}{see \code{\link{par}} function.}
  \item{vectors}{Whether vectors representing direction of photoreceptor outputs should be plotted.}
  \item{...}{Other arguments passed to \code{\link{plot}}.}   
}
\details{
The original model is available for tetrachromatic animals only. Trichromatic version is implemented in \code{colourvision} based on Pike (2012) formula. The triangle is built either with a vector length = 0.75 or a edge length = sqrt(3/2), to match the tetrahedron proposed by Endler and Mielke (2005). Doris Gomez derived a trichromatic version which is available in software AVICOL (Gomez, 2006) and was previously implemented here (\code{colouvision v0.1}).
}
\references{
Endler, J. A., and P. Mielke. 2005. Comparing entire colour patterns as birds see them. Biol J Linn Soc 86:405-431.

Pike, T.W. 2012. Generalised chromaticity diagrams for animals with n-chromatic colour vision. Journal of Insect Behavior 255: 277-286.

Gomez, D. 2006. AVICOL, a program to analyse spectrometric data. Last update october 2011. Free executable available at:\cr
\url{http://sites.google.com/site/avicolprogram/}  or from the author at \email{dodogomez@yahoo.fr}
}
\author{
Felipe M. Gawryszewski \email{f.gawry@gmail.com}
}
\seealso{
\code{\link{EMtetrahedron}},  \code{\link{EMmodel}}
}
\examples{
EMtriangle(x=0,y=0, type="length", pch=16, col="red")
}