\name{RNLthres}
\alias{RNLthres}
\title{
Colour thresholds based on the Receptor Noise Limited Model (Vorobyev & Osorio 1998).
}
\description{
Colour thresholds based on receptor noise for any number of photoreceptor types (Vorobyev & Osorio 1998).
}
\usage{
RNLthres(photo=ncol(C)-1, Rb, I, C, noise=TRUE, v=NA, n=NA, e=NA,
         interpolate=TRUE, nm=seq(300,700,1))
}
\arguments{
  \item{photo}{
Number of photoreceptor types. Model accepts any number of photoreceptor types (\code{>=2}). For instance, dichromatic: \code{photo=2}; trichromatic: \code{photo=3}; tetrachromatic: \code{photo=4}, etc. Default gets number of photoreceptor types from \code{C} argument.
}
  \item{Rb}{
Background reflectance. A data frame with two columns only: first column corresponding to wavelength values and second column with reflectance values.
}
  \item{I}{
Irradiance spectrum. A data frame with two columns only: first column corresponding to wavelength values and second column with irradiance values. Irradiance values must be in quantum flux units.
}
  \item{C}{
Photoreceptor sensitivity curves, from lowest to longest lambda-max. A data frame: first column corresponding to wavelength values and following columns with photoreceptor sensitivity values (see function \code{\link{photor}}). 
}
  \item{noise}{
Logical. Whether receptor noise is provided (\code{noise = TRUE}) or calculated from photoreceptor relative abundances (\code{noise = FALSE}).
}
  \item{e}{
Receptor noise of each photoreceptor type. Used when \code{noise = TRUE}
}
  \item{n}{
Relative number of each photoreceptor type in the retina. Usually increases with lambda-max. Used to calculate \code{e} when \code{noise = FALSE}.
}
  \item{v}{
Noise-to-signal ratio of a single photoreceptor. Used to calculate \code{e} when \code{noise = FALSE}.
}
  \item{interpolate}{
Whether data files should be interpolated before further calculations. See \code{\link{approx}}.}
  \item{nm}{
A sequence of numeric values specifying where interpolation is to take place. See \code{\link{approx}}.}
}
\value{
A data.frame with the following columns:
\item{nm}{Wavelength in nm.}
\item{T}{Colour threshold value.}
\item{S}{Log of sensitivity value (inverse of threshold).}
}

\details{
Colour thresholds based on receptor noise limited model as in Vorobyev and Osorio (1998). In lack of a direct measurement, receptor noise (\eqn{e_i}{ei}) can be estimated by the relative abundance of photoreceptor types in the retina, and a measurement of a single photoreceptor noise-to-signal ratio:
\deqn{e_i=\frac{\nu}{\sqrt{\eta _i}}}{ei=v/sqrt(n)}
where \eqn{\nu}{v} is the noise-to-signal ratio of a single photoreceptor, and \eqn{\eta}{n} is the relative abundance of photoreceptor i in the retina. Alternatively, noise may be dependent of the intensity, but this possibility is not implement in \code{colourvision} yet. Noise dependent of intensity usually holds for low light conditions only (Vorobyev et al. 1998). 
}
\references{
Vorobyev, M., and D. Osorio. 1998. Receptor noise as a determinant of colour thresholds. Proceedings of the Royal Society B 265:351-358.
}

\author{
Felipe M. Gawryszewski \email{f.gawry@gmail.com}
}
\seealso{
\code{\link{photor}}, \code{\link{RNLmodel}}
}
\examples{
###Bee photoreceptors normalized to max=1.
data("bee")
C<-bee
C[,2]<-C[,2]/max(C[,2])
C[,3]<-C[,3]/max(C[,3])
C[,4]<-C[,4]/max(C[,4])

##Grey background:
Rb <- data.frame(300:700, rep(0.1, length(300:700)))

## CIE D65 illuminant:
data("D65")

#Thresholds
thres<-RNLthres(photo=3, Rb=Rb, I=D65, C=C,
       noise=TRUE, e = c(0.13, 0.06, 0.12))

plot(thres)     
}