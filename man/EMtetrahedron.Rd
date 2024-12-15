\name{EMtetrahedron}
\alias{EMtetrahedron}
\title{
Endler and Mielke (2005) tetrahedron colour space}
\description{
Plots Endler and Mielke (2005) tetrahedron colour space for tetrachromatic animals.}
\usage{
EMtetrahedron(x, y, z, s.col = "red", f.col = "black",
              vnames = c("u","s","m","l"), type = "p",
              radius = 0.01, add = F, xlab = "",
              ylab = "", zlab = "", box = F, axes = F,
              ylim = c(-0.75, 0.75), xlim = c(-0.75, 0.75),
              zlim = c(-0.75, 0.75), aspect = T, vectors=FALSE, ...)
}
\arguments{
  \item{x}{x coordinate of points to be plotted}
  \item{y}{y coordinate of points to be plotted}
  \item{z}{z coordinate of points to be plotted}
  \item{s.col}{Colour to be used for plotted items. See plot3d in the rgl package.}
  \item{f.col}{Colour of tetrahedron lines. See plot3d in the rgl package.}
  \item{vnames}{Vector names.}
  \item{type}{see plot3d function in the rgl package.}
  \item{radius}{see plot3d function in the rgl package.}
  \item{add}{see plot3d function in the rgl package.}
  \item{xlab}{see plot3d function in the rgl package.}
  \item{ylab}{see plot3d function in the rgl package.}
  \item{zlab}{see plot3d function in the rgl package.}
  \item{box}{see plot3d function in the rgl package.}
  \item{axes}{see plot3d function in the rgl package.}
  \item{ylim}{see plot3d function in the rgl package.}
  \item{xlim}{see plot3d function in the rgl package.}
  \item{zlim}{see plot3d function in the rgl package.}
  \item{aspect}{see plot3d function in the rgl package.}
  \item{vectors}{Whether vectors representing direction of photoreceptor outputs should be plotted.}
  \item{...}{Other arguments passed to function plot3d in the rgl package.}
}
\references{
Endler, J. A., and P. Mielke. 2005. Comparing entire colour patterns as birds see them. Biol J Linn Soc 86:405-431.
}
\author{
Felipe M. Gawryszewski \email{f.gawry@gmail.com}
}
\seealso{
\code{\link{EMtriangle}}, \code{\link{EMmodel}}
}