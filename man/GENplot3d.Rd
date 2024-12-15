\name{GENplot3d}
\alias{GENplot3d}
\title{
Generic model colour space 3D plot
}
\description{
Plots models based on the GENmodel( ) function for tetrachromatic animals.
}
\usage{
GENplot3d(model, col.names=c("X1","X2","X3"),
          vectors=TRUE, vnames=TRUE, vsize="auto",
          xlab="x", ylab="y", zlab="z",
          xlim="auto", ylim="auto", zlim="auto", asp=1, ...)
}
\arguments{
  \item{model}{Output of a colour vision model.}
  \item{col.names}{Column names to be plotted.}
  \item{vectors}{Whether vectors representing direction of photoreceptor outputs should be plotted.}
  \item{vnames}{Whether vector names should be plotted.}
  \item{vsize}{Length of vectors to be plotted. Default calculates length automatically.}
  \item{xlab}{see par3d function in the rgl package.}
  \item{ylab}{see par3d function in the rgl package.}
  \item{zlab}{see par3d function in the rgl package.}
  \item{xlim}{x-axis range. Default calculates range automatically. See par3d function in the rgl package.}
  \item{ylim}{y-axis range. Default calculates range automatically. See par3d function in the rgl package.}
  \item{zlim}{z-axis range. Default calculates range automatically. See par3d function in the rgl package.}
  \item{asp}{see plot3d function in the rgl package.}
  \item{...}{Other arguments passed to function plot3d in the rgl package.}
}
\author{
Felipe M. Gawryszewski \email{f.gawry@gmail.com}
}
\seealso{
\code{\link{CTTKhexagon}}, \code{\link{CTTKhexagon3D}}, \code{\link{EMtriangle}}, \code{\link{EMtetrahedron}}, \code{\link{RNLplot}}, \code{\link{RNLplot3d}}, \code{\link{GENplot}}, \code{\link{plot.colourvision}}, \code{\link{plot3d.colourvision}}
}