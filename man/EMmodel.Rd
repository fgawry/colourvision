\name{EMmodel}
\alias{EMmodel}
\title{
Endler and Mielke (2005) colour vision model
}
\description{
Endler and Mielke (2005) colour vision model extended to animals with any number of photoreceptor types.
}
\usage{
EMmodel(photo = ncol(C)-1, type="length", R, I, Rb, C,
        interpolate=TRUE, nm=seq(300,700,1))
}
\arguments{
  \item{photo}{
Number of photoreceptor types. Model accepts any number of photoreceptor types (\code{>=2}). For instance, dichromatic: \code{photo=2}; trichromatic: \code{photo=3}; tetrachromatic: \code{photo=4}, etc. Default gets number of photoreceptor types from \code{C} argument.}
  \item{type}{Whether the colour space should be built with a fixed vector length (\code{type="length"}), or a fixed edge length (distance between vertices; \code{type="edge")}.
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
The original model is available for tetrachromatic animals only. In \code{colourvision}, the model was extended to any number of photoreceptors types (see also Pike 2012 formula).

First, relative quantum catches are log-transformed:
\deqn{f_{i} = \ln{q_i}}{fi = ln(qi)}
where \eqn{q_i}{qi} is the relative quantum catch of photoreceptor type i, given by \code{\link{Qr}}. The model uses only relative output values, so that photoreceptor outputs are given by:
\deqn{E_i = \frac{f_i}{\sum_{i=1}^n f_i}}{Ei = fi/sum(f1+f2+f3+...fn)}

For tetrachromatic vision (Endler and Mielke 2005):
\deqn{X1 = \sqrt{\frac{3}{2}}(\frac{1-2E_2-E_3-E_1}{2})}{X1 = sqrt(3/2)*((1-2*E2-E3-E1)/2)}
\deqn{X2 = \frac{-1+3E_3+E_1}{2\sqrt{2}}}{X2 = (-1+3*E3+E1)/2*sqrt(2)}
\deqn{X3 = E_1-\frac{1}{4}}{X3 = E1-0.25}

Tetrachromatic chromaticity diagram (tetrahedron) in Endler and Mielke (2005) has a vector of length = 0.75 and and edge length = sqrt(3/2). The chromaticity coordinates for other colour spaces may preserve either the same vector length or edge length. 

For instance, for dichromatic vision coordinate (X1) in the colour space preserving the same vector length is found by:
\deqn{X1 = \frac{3}{4}(E_2-E_1)}{X1 = (3/4)*(E2-E1)}
Whereas for trichromatic vision coordinates (X1 and X2) are found by:
\deqn{X1 = \frac{3\sqrt{3}}{8}(E_2-E_1)}{X1 = (3*sqrt(3)/8)*(E2-E1)}
\deqn{X2 = \frac{3}{4}(E_3-\frac{E_2+E_1}{2})}{X2 = (3/4)*(E3-(E2+E1)/2)}
}

\value{
\item{Qri}{Photoreceptor photon catch values after the von Kries transformation (see function \code{\link{Qr}}).}
\item{Ei}{Photoreceptor outputs after conversion to relative values.}
\item{Xi}{Coordinates in the colour space.}
\item{deltaS}{Euclidean distance to the origin of the colour space. It represents the conspicuousness of the stimulus (\code{R}) in relation to the background (\code{Rb}).}
}
\references{
Endler, J. A., and P. Mielke. 2005. Comparing entire colour patterns as birds see them. Biol J Linn Soc 86:405-431.

Pike, T.W. 2012. Generalised chromaticity diagrams for animals with n-chromatic colour vision. Journal of Insect Behavior 255: 277-286.
}
\author{
Felipe M. Gawryszewski \email{f.gawry@gmail.com}
}
\seealso{
\code{\link{EMline}}, \code{\link{EMtriangle}}, \code{\link{EMtetrahedron}}, \code{\link{photor}}, \code{\link{CTTKmodel}}, \code{\link{RNLmodel}}, \code{\link{GENmodel}}
}
\examples{
##Photoreceptor sensitivity curves
##with lambda max at 350nm, 450nm and 550nm:
C<-photor(lambda.max=c(350,450,550))

##Gray background
##with 7 percent reflectance from 300 to 700nm:
Rb <- data.frame(300:700, rep(7, length(300:700)))

## Read CIE D65 standard illuminant
data("D65")

##Reflectance data
## with a sigmoid spectrum and midpoint at 500nm and 550 nm
R1<-logistic(x=seq(300,700,1), x0=500, L=50, k=0.04)
R2<-logistic(x=seq(300,700,1), x0=550, L=50, k=0.04)
R<-cbind(R1, R2[,2])
R[,2]<-R[,2]+10
R[,3]<-R[,3]+10

## Run model 
model<-EMmodel(photo=3, type="edge",
       R=R, I=D65, Rb=Rb, C=C)

plot(model)
}