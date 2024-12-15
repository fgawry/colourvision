\name{plot3d.colourvision}
\alias{plot3d.colourvision}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot colour vision models into 3D chromaticity diagrams.
}
\description{
'plot3d' method for objects of class \code{colourvision}. Plotting method for animals with four photoreceptor types.
}
\usage{
\method{plot3d}{colourvision}(x, ...)
}

\arguments{
  \item{x}{
Object of class "\code{colourvision}".
}
  \item{\dots}{
Additional arguments passed to the plot3d function. See \code{\link{CTTKhexagon3D}} for Chittka (1992) model (\code{\link{CTTKmodel}}) plotting arguments; \code{\link{EMtetrahedron}} for Endler and Mielke (2005) model (\code{\link{EMmodel}}) plotting arguments;  \code{\link{RNLplot3d}} for RNL plotting arguments; and \code{\link{GENplot3d}} for generic model (\code{\link{GENmodel}}) plotting arguments.
}
}
\references{
Chittka, L. 1992. The colour hexagon: a chromaticity diagram based on photoreceptor excitations as a generalized representation of colour opponency. J Comp Physiol A 170:533-543.

Endler, J. A., and P. Mielke. 2005. Comparing entire colour patterns as birds see them. Biol J Linn Soc 86:405-431.

Thery, M., and J. Casas. 2002. Predator and prey views of spider camouflage. Nature 415:133-133.
}

\author{
Felipe M. Gawryszewski \email{f.gawry@gmail.com}
}
\seealso{
\code{\link{plot.colourvision}}, \code{\link{EMtetrahedron}}, \code{\link{CTTKhexagon3D}}, \code{\link{EMmodel}}, \code{\link{CTTKmodel}}, \code{\link{RNLmodel}}
}