\name{RNLachrom}
\alias{RNLachrom}
\title{
Weber achromatic contrast for the Receptor Noise Limited Model (Vorobyev & Osorio 1998)
}
\description{
Weber achromatic contrast for the Receptor noise limited model (Vorobyev & Osorio 1998; Vorobyev et al. 1998).
}
\usage{
RNLachrom(R1, R2=Rb, Rb, I, C, e,
         interpolate = TRUE, nm = seq(300, 700, 1))
}
\arguments{
\item{R1}{
Reflectance of observed objects. A data frame with first column corresponding to wavelength values and following columns with reflectance values. \code{R1} must be in the same scale as \code{R2} and \code{Rb} (percentage or proportion).
}
\item{R2}{
Reflectance to be compared against R1. \code{R2=Rb} if comparison is to be made against the background. A data frame with two columns only: first column corresponding to wavelength values and second column with reflectance values. \code{R2} must be in the same scale as \code{R1} and \code{Rb} (percentage or proportion).
}
\item{Rb}{
Background reflectance. A data frame with two columns only: first column corresponding to wavelength values and second column with reflectance values. \code{Rb} must be in the same scale as \code{R1} and \code{R2} (percentage or proportion).
}
\item{I}{
Irradiance spectrum. A data frame with two columns only: first column corresponding to wavelength values and second column with irradiance values. Irradiance values must be in quantum flux units.
}
\item{C}{
Photoreceptor sensitivity curves, from lowest to longest lambda-max. A data frame: first column corresponding to wavelength values and following columns with photoreceptor sensitivity values (see function \code{\link{photor}}). 
}
\item{e}{
Noise of the photoreceptor.
}
\item{interpolate}{
Whether data files should be interpolated before further calculations. See \code{\link{approx}}.}
\item{nm}{
A sequence of numeric values specifying where interpolation is to take place. See \code{\link{approx}}.
}

}
\value{
\item{e1}{Photoreceptor noise.}
\item{Qr1_R1}{Photoreceptor photon catch value from \code{R1}, after the von Kries transformation (see function \code{\link{Qr}}).}
\item{Qr1_R2}{Photoreceptor photon catch value from \code{R2}, after the von Kries transformation (see function \code{\link{Qr}}).}
\item{E1_R1}{Photoreceptor output from the stimulus (\code{R1})}
\item{E1_R2}{Photoreceptor output from R2}
\item{deltaS}{Euclidean distance from R1 to R2. It represents the achromatic contrast of the stimulus (\code{R1}) in relation to the background when \code{R1=Rb}.}
}

\details{

The Weber achromatic contrast for a single photoreceptor is calculated by:
\deqn{\Delta S = |\frac{\ln(Qr_1)-\ln(Qr_2)}{e}|}{}

where \eqn{Qr_1} and \eqn{Qr_2} are the relative photoreceptor quantum catches from stimulus 1 (\code{R1}) and stimulus 2 (\code{R2}).

Noise may be dependent of the intensity, but this possibility is not implement in \code{colourvision} yet. Noise dependent of intensity usually holds for low light conditions only (Vorobyev et al. 1998). 
}

\references{
Vorobyev, M., and D. Osorio. 1998. Receptor noise as a determinant of colour thresholds. Proceedings of the Royal Society B 265:351-358.

Vorobyev, M., D. Osorio, A. T. D. Bennett, N. J. Marshall, and I. C. Cuthill. 1998. Tetrachromacy, oil droplets and bird plumage colours. J Comp Physiol A 183:621-633.
}

\author{
Felipe M. Gawryszewski \email{f.gawry@gmail.com}
}
\seealso{
\code{\link{RNLmodel}}, \code{\link{photor}}, \code{\link{RNLthres}}, \code{\link{CTTKmodel}}, \code{\link{EMmodel}}, \code{\link{GENmodel}}
}
\examples{
#1
## Photoreceptor sensitivity spectra
##with lambda max at 350nm, 450nm and 550nm:
C<-photor(lambda.max=c(350))

##Grey background
##with 7 percent reflectance from 300 to 700nm:
Rb <- data.frame(300:700, rep(7, length(300:700)))

## Read CIE D65 standard illuminant:
data("D65")

##Reflectance data of R1 and R2
R1.1<-logistic(x=seq(300,700,1), x0=500, L=50, k=0.04)
R1.2<-logistic(x=seq(300,700,1), x0=400, L=50, k=0.04)
w<-R1.1[,1]
R1.1<-R1.1[,2]+10
R1.2<-R1.2[,2]+10
R1<-data.frame(w=w, R1.1=R1.1, R1.2=R1.2)

R2<-logistic(x=seq(300,700,1), x0=550, L=50, k=0.04)
R2[,2]<-R2[,2]+10

plot(R1[,c(1,2)],type="l",
     ylim=c(0,60))
lines(R1[,c(1,3)])
lines(R2[,c(1,2)],col="red")
lines(Rb,col="green")

## Run model 
RNLachrom(R1=R1, R2=R2, Rb=Rb, I=D65, C=C,
          e = 0.16)

}