\name{Qr}
\alias{Qr}
\title{
Photoreceptor relative quantum catch}
\description{
von Kries transformation. Photoreceptors are assumed to be adapted to the background. This function is used internally in colour vision models.
}
\usage{
Qr(R, I, Rb, C, interpolate, nm)
}
\arguments{
  \item{R}{
Reflectance of observed object. A data frame with two columns only: first column corresponding to wavelength values and second column with reflectance values.
}
  \item{I}{
Irradiance spectrum. A data frame with two columns only: first column corresponding to wavelength values and second column with irradiance values. Irradiance values must be in quantum flux units.
}
  \item{Rb}{
Background reflectance. A data frame with two columns only: first column corresponding to wavelength values and second column with reflectance values. Photoreceptors are assumed to be adapted to the background reflectance.
}
\item{C}{
Photoreceptor sensitivity curve. A data frame with two columns only: first column corresponding to wavelength values and second column with photoreceptor absorbance values.
}
  \item{interpolate}{
Whether data files should be interpolated before further calculations. See \code{\link{approx}}.}
  \item{nm}{
A sequence of numeric values specifying where interpolation is to take place. See \code{\link{approx}}.}
}
\details{For the von Kries transformation, first the quantum catches of the observed reflectance and the environmental background are calculated (see \code{\link{Q}}). Then:
\deqn{qi = \frac{Q_i}{Q_{bi}}}{qi = Qi/Qbi}
where \eqn{Q_i}{Qi} is the quantum catch arising from the observed object and \eqn{Q_{bi}}{Qbi} is the quantum catch from the background, for each one of the photoreceptor types (i).
}
\value{
Photoreceptor relative quantum catch.
}
\references{
Backhaus, W. 1991. Color opponent coding in the visual system of the honeybee. Vision Res 31:1381-1397.

Chittka, L. 1992. The colour hexagon: a chromaticity diagram based on photoreceptor excitations as a generalized representation of colour opponency. J Comp Physiol A 170:533-543.

Endler, J. A., and P. Mielke. 2005. Comparing entire colour patterns as birds see them. Biol J Linn Soc 86:405-431.

Vorobyev, M., and D. Osorio. 1998. Receptor noise as a determinant of colour thresholds. Proceedings of the Royal Society B 265:351-358.
}
\author{
Felipe M. Gawryszewski \email{f.gawry@gmail.com}
}
\seealso{
\code{\link{CTTKmodel}}, \code{\link{EMmodel}}, \code{\link{RNLmodel}}, \code{\link{RNLthres}}, \code{\link{GENmodel}}
}