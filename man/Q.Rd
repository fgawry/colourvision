\name{Q}
\alias{Q}
\title{
Total photon capture}
\description{
Total photoreceptor photon capture for a given irradiance, reflectance and photoreceptor sensitivity curve.  This function is used internally in colour vision models.}
\usage{
Q(R,I,C,interpolate,nm)
}
\arguments{
  \item{R}{
Reflectance of observed object. A data frame with two columns only: first column corresponding to wavelength values and second column with reflectance values.
}
  \item{I}{
Irradiance spectrum. A data frame with two columns only: first column corresponding to wavelength values and second column with irradiance values. Irradiance values must be in quantum flux units.
}
  \item{C}{
Photoreceptor sensitivity curve. A data frame with two columns only: first column corresponding to wavelength values and second column with photoreceptor absorbance values.
}
  \item{interpolate}{
Whether data files should be interpolated before further calculations. See \code{\link{approx}}.}
  \item{nm}{
A sequence of numeric values specifying where interpolation is to take place. See \code{\link{approx}}.}
}
\value{
Gives the total photoreceptor photon capture.
}
\references{
Backhaus, W., and R. Menzel. 1987. Color distance derived from a receptor model of color vision in the honeybee. Biological Cybernetics 55:321-331.

Chittka, L. 1992. The colour hexagon: a chromaticity diagram based on photoreceptor excitations as a generalized representation of colour opponency. J Comp Physiol A 170:533-543.

Endler, J. A., and P. Mielke. 2005. Comparing entire colour patterns as birds see them. Biol J Linn Soc 86:405-431.

Vorobyev, M., and D. Osorio. 1998. Receptor noise as a determinant of colour thresholds. Proceedings of the Royal Society B 265:351-358.
}
\author{
Felipe M. Gawryszewski \email{f.gawry@gmail.com}
}
\seealso{
\code{\link{Qr}}, \code{\link{CTTKmodel}}, \code{\link{EMmodel}}, \code{\link{RNLmodel}}, \code{\link{RNLthres}},\code{\link{GENmodel}}
}