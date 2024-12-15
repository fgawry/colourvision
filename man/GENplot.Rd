\name{GENplot}
\alias{GENplot}
\title{
Generic model colour space 2D and 1D plot}
\description{
Plots models based on the GENmodel( ) function for trichromatic and dichromatic animals.
}
\usage{
GENplot(model, photo, col.names=c("X1","X2"),
        vectors=TRUE, vnames=TRUE, vsize="auto",
        ylab="y", xlab="x", xlim="auto", ylim="auto", asp=1, ...)
}
\arguments{
  \item{model}{Output of a colour vision model.}
  \item{photo}{Number of photoreceptor types.}
  \item{col.names}{Column names to be plotted.}
  \item{vectors}{Whether vectors representing direction of photoreceptor outputs should be plotted.}
  \item{vnames}{Whether vector names should be plotted.}
  \item{vsize}{Length of vectors to be plotted. Default calculates length automatically.}
  \item{ylab}{y-axis range. Default calculates range automatically. See \code{\link{par}} function.}
  \item{xlab}{x-axis range. Default calculates range automatically. See \code{\link{par}} function.}
  \item{xlim}{see \code{\link{par}} function.}
  \item{ylim}{see \code{\link{par}} function.}
  \item{asp}{see \code{\link{plot}} function.}
  \item{...}{Other arguments passed to \code{\link{plot}} function.}
}
\author{
Felipe M. Gawryszewski \email{f.gawry@gmail.com}
}
\seealso{
\code{\link{CTTKhexagon}}, \code{\link{CTTKhexagon3D}}, \code{\link{EMtriangle}}, \code{\link{EMtetrahedron}}, \code{\link{RNLplot}}, \code{\link{RNLplot3d}}, \code{\link{GENplot3d}}, \code{\link{plot.colourvision}}, \code{\link{plot3d.colourvision}}
}