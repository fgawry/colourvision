\name{CTTKmodel}
\alias{CTTKmodel}
\title{
Chittka (1992) colour vision model
}
\description{
Chittka (1992) colour hexagon extended to animals with any number of photoreceptors types.
}
\usage{
CTTKmodel(photo=ncol(C)-1, R, I, Rb, C,
          interpolate=TRUE, nm=seq(300,700,1))
}
\arguments{
  \item{photo}{
Number of photoreceptor types. Model accepts any number of photoreceptor types (\code{>=2}). For instance, dichromatic: \code{photo=2}; trichromatic: \code{photo=3}; tetrachromatic: \code{photo=4}, etc. Default gets number of photoreceptor types from \code{C} argument.
}
  \item{R}{
Reflectance of observed objects. A data frame with first column corresponding to wavelength values and following columns with reflectance values. \code{R} must be in the same scale as \code{Rb} (percentage or proportion).
}
  \item{I}{
Irradiance spectrum. A data frame with two columns only: first column corresponding to wavelength values and second column with irradiance values. Irradiance values must be in quantum flux units.
}
  \item{Rb}{
Background reflectance. A data frame with two columns only: first column corresponding to wavelength values and second column with reflectance values. \code{Rb} must be in the same scale as \code{R} (percentage or proportion).
}
  \item{C}{
Photoreceptor sensitivity curves, from lowest to longest lambda-max. A data frame: first column corresponding to wavelength values and following columns with photoreceptor sensitivity values (see function \code{\link{photor}}). 
}
  \item{interpolate}{
Whether data files should be interpolated before further calculations. See \code{\link{approx}}.}
  \item{nm}{
A sequence of numeric values specifying where interpolation is to take place. See \code{\link{approx}}.}

}
\details{
The original model is available for trichromatic animals only. Thery and Casas (2002) derived a version for tetrachromatic animals which is implemented here. In \code{colourvision}, this model was extended to any number of photoreceptors types (Gawryszewski 2018; see also Pike 2012). The colour hexagon in Chittka (1992) has a vector of length = 1.0 The chromaticity coordinates in \code{colourvision} preserve the same vector length.

Photoreceptor outputs (\eqn{E_i}{Ei}) are calculated by:
\deqn{E_i = \frac{q_i}{q_i+1}}{Ei = qi/(qi+1)}
where \eqn{q_i}{qi} is given by \code{\link{Qr}}.

Then, for trichromatic vision, coordinates in the colour space are found by (Chittka 1992):
\deqn{X_1 = \frac{\sqrt{3}}{2}(E_3-E_1)}{X1 = (sqrt(3)/2)*(E3-E1)}
\deqn{X_2 = E_2-\frac{1}{2}(E_1+E_3)}{X2 = E2-0.5*(E1+E3)}
For tetrachromatic vision (Thery and Casas 2002):
\deqn{X_1 = \frac{\sqrt{3}\sqrt{2}}{3}(E_3-E_4)}{X1 = ((sqrt(3)*sqrt(2))/3)*(E3-E4)}
\deqn{X_2 = E_1-\frac{1}{3}(E_2+E_3+E_4)}{X2 = E1-(1/3)*(E2+E3+E4)}
\deqn{X_3 = \frac{2\sqrt{2}}{3}(\frac{1}{2}(E_3+E_4)-E_2)}{X3 = (2*sqrt(2)/3)*(0.5*(E3+E4)-E2)}
For a pentachromatic animal following the same vector length:
\deqn{X_1 = \frac{5}{2\sqrt{2}\sqrt{5}}(E_2-E_1)}{X1 = (5/(2*sqrt(2)*sqrt(5)))*(E2-E1)}
\deqn{X_2 = \frac{5\sqrt{2}}{2\sqrt{3}\sqrt{5}}(E_3-\frac{E_1+E_2}{2})}{X2 = (5*sqrt(2)/(2*sqrt(3)*sqrt(5)))*(E3-(E1+E2)/2)}
\deqn{X_3 = \frac{5\sqrt{3}}{4\sqrt{5}}(E_4-\frac{E_1+E_2+E_3}{3})}{X3 = (5*sqrt(3)/(4*sqrt(5)))*(E4-(E1+E2+E3)/3)}
\deqn{X_4 = E_5-\frac{E1+E2+E3+E4}{4}}{X4 = E5-(E1+E2+E3+E4)/4}

}
\value{
\item{Qri}{Photoreceptor photon catch values after the von Kries transformation (see function \code{\link{Qr}}).}
\item{Ei}{Photoreceptor output values. Values can vary from 0 to 1.}
\item{Xi}{Coordinates in the colour space.}
\item{deltaS}{Euclidean distance to the origin of the colour space. It represents the conspicuousness of the stimulus (\code{R}) in relation to the background (\code{Rb}).}
}
\references{
Chittka, L. 1992. The colour hexagon: a chromaticity diagram based on photoreceptor excitations as a generalized representation of colour opponency. J Comp Physiol A 170:533-543.

Gawryszewski, F.M. 2018. Colour vision models: Some simulations, a general n-dimensional model, and the colourvision R package. Ecology and Evolution, 10.1002/ece3.4288.

Pike, T.W. 2012. Generalised chromaticity diagrams for animals with n-chromatic colour vision. Journal of Insect Behavior 255: 277-286.

Thery, M., and J. Casas. 2002. Predator and prey views of spider camouflage. Nature 415:133-133.
}
\author{
Felipe M. Gawryszewski \email{f.gawry@gmail.com}
}
\seealso{
\code{\link{CTTKhexagon}}, \code{\link{CTTKhexagon3D}}, \code{\link{photor}}, \code{\link{RNLmodel}}, \code{\link{EMmodel}}, \code{\link{deltaS}}
}
\examples{
##Photoreceptor sensitivity curves
##with lambda max at 350nm, 450nm and 550nm:
C<-photor(lambda.max=c(350,450,550))

## Grey background
## with 10 percent reflectance from 300 to 700nm:
Rb <- data.frame(300:700, rep(10, length(300:700)))

## Read CIE D65 standard illuminant
data("D65")

## Reflectance data
## with a sigmoid spectrum and midpoint at 500nm
R<-logistic(x=seq(300,700,1), x0=500, L=50, k=0.04)

## Run model 
model<-CTTKmodel(photo=3, R=R, I=D65,
    Rb=Rb, C=C)
 
#plot   
plot(model)
}