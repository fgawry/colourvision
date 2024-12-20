\name{RNLmodel}
\alias{RNLmodel}
\title{
Receptor Noise Limited Models (Vorobyev & Osorio 1998)
}
\description{
Receptor noise limited colour vision models (Vorobyev & Osorio 1998; Vorobyev et al. 1998) extended to any number of photoreceptor types.
}
\usage{
RNLmodel(model = c("linear", "log"), photo=ncol(C)-1,
         R1, R2=Rb, Rb, I, C,
         noise = FALSE, v=NA, n=NA, e=NA,
         interpolate = TRUE, nm = seq(300, 700, 1),
         coord="colourvision")
}
\arguments{
\item{model}{Linear (\code{"linear"}; Vorobyev & Osorio, 1998), or log-linear  (\code{"log"}; Vorobyev et al. 1998) RNLmodel version. \code{model="log"} is preferred.

\code{model="linear"} can be used only if comparison is to be made between two very similar colours.
}
  \item{photo}{
Number of photoreceptor types. Model accepts any number of photoreceptor types (\code{>=2}). For instance, dichromatic: \code{photo=2}; trichromatic: \code{photo=3}; tetrachromatic: \code{photo=4}, etc. Default gets number of photoreceptor types from \code{C} argument. If only a single photoreceptor is provided (e.g., \code{photo = 1}, uses the \code{\link{RNLachrom}} function for achromatic contrast calculations.
}
  \item{R1}{
Reflectance of observed objects. A data frame with first column corresponding to wavelength values and following columns with reflectance values. \code{R1} must be in the same scale as \code{R2} and \code{Rb}(percentage or proportion).
}
  \item{R2}{
Reflectance to be compared against R1. \code{R2=Rb} if comparison is to be made against the background. A data frame with two columns only: first column corresponding to wavelength values and second column with reflectance values. \code{R2} must be in the same scale as \code{R1} and \code{Rb}(percentage or proportion).
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
Noise-to-signal ratio of a single photoreceptor. Either a single value or a vector of same length as \code{n}. If the former, assumes that the \code{v} value refers to the most common value in \code{n}. If the latter, the vector must have only one numeric value, and the remaining as NA (e.g., \code{v = c(NA, 0.1, NA)}).
}
  \item{interpolate}{
Whether data files should be interpolated before further calculations. See \code{\link{approx}}.}
  \item{nm}{
A sequence of numeric values specifying where interpolation is to take place. See \code{\link{approx}}.}
  \item{coord}{Whether colour locous coordinates should be calculated by the method (\code{coord="colourvision"}) by Gawryszewski (2018), or by alternative (\code{coord="alternative"}; available for 2-4 photoreceptor types only) methods (Hempel de Ibarra et al. 2001; Renoult et al. 2017).}

}
\value{
\item{ei}{Noise of photoreceptor channels.}
\item{Qri_R1}{Photoreceptor photon catch values from \code{R1}, after the von Kries transformation (see function \code{\link{Qr}}).}
\item{Qri_R2}{Photoreceptor photon catch values from \code{R2}, after the von Kries transformation (see function \code{\link{Qr}}).}
\item{Ei_R1}{Photoreceptor outputs from the stimulus (\code{R1})}
\item{Ei_R2}{Photoreceptor outputs from R2}
\item{Xi_R1}{Coordinates in the colour space for R1}
\item{Xi_R2}{Coordinates in the colour space for R2. Equals zero when \code{R1=Rb}}
\item{deltaS}{Euclidean distance from R1 to R2. It represents the conspicuousness of the stimulus (\code{R1}) in relation to the background when \code{R1=Rb}.}
}

\details{
The receptor noise limited model was originally developed to calculate \eqn{\Delta S}{deltaS} between two reflectance curves directly, without finding colour locus coordinates (e.g. \code{x,y}; Vorobyev and Osorio 1998). This function uses later formulae to find colour loci in a chromaticity diagram (similarly to Hempel de Ibarra et al. 2001; Renoult et al. 2015).

In lack of a direct measurement, receptor noise (\eqn{e_i}{ei}) can be estimated by the relative abundance of photoreceptor types in the retina, and a measurement of a single photoreceptor noise-to-signal ratio:
\deqn{e_i=\frac{\nu}{\sqrt{\eta _i}}}{ei=v/sqrt(n)}
where \eqn{\nu}{v} is the noise-to-signal ratio of a single photoreceptor, and \eqn{\eta}{n} is the relative abundance of photoreceptor i in the retina. Alternatively, noise may be dependent of the intensity, but this possibility is not implement in \code{colourvision} yet. Noise dependent of intensity usually holds for low light conditions only (Vorobyev et al. 1998). 
}
\references{
Gawryszewski, F.M. 2018. Colour vision models: Some simulations, a general n-dimensional model, and the colourvision R package. Ecology and Evolution, 10.1002/ece3.4288.

Hempel de Ibarra, N., M. Giurfa, and M. Vorobyev. 2001. Detection of coloured patterns by honeybees through chromatic and achromatic cues. J Comp Physiol A 187:215-224.

Renoult, J. P., A. Kelber, and H. M. Schaefer. 2017. Colour spaces in ecology and evolutionary biology. Biol Rev Camb Philos Soc, doi: 10.1111/brv.12230

Vorobyev, M., and D. Osorio. 1998. Receptor noise as a determinant of colour thresholds. Proceedings of the Royal Society B 265:351-358.

Vorobyev, M., D. Osorio, A. T. D. Bennett, N. J. Marshall, and I. C. Cuthill. 1998. Tetrachromacy, oil droplets and bird plumage colours. J Comp Physiol A 183:621-633.
}

\author{
Felipe M. Gawryszewski \email{f.gawry@gmail.com}
}
\seealso{
\code{\link{RNLachrom}}, \code{\link{photor}}, \code{\link{RNLthres}}, \code{\link{CTTKmodel}}, \code{\link{EMmodel}}, \code{\link{GENmodel}}
}
\examples{
#1
## Photoreceptor sensitivity spectra
##with lambda max at 350nm, 450nm and 550nm:
C<-photor(lambda.max=c(350,450,550))

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

## Run model 
model<-RNLmodel(photo=3, model="log",
       R1=R1, R2=R2, Rb=Rb, I=D65, C=C,
       noise=TRUE, e = c(0.13, 0.06, 0.12))
       
#plot       
plot(model)

#2
#Pentachromatic animal
## Photoreceptor sensitivity spectra
##with lambda max at 350,400,450,500,and 550nm:
C<-photor(lambda.max=c(350,400,450,500,550))

##Grey background
##with 7 percent reflectance from 300 to 700nm:
Rb <- data.frame(300:700, rep(7, length(300:700)))

## Read CIE D65 standard illuminant:
data("D65")

##Reflectance data of R1
R1<-logistic(x=seq(300,700,1), x0=500, L=50, k=0.04)
R1[,2]<-R1[,2]+10

#RNL model
RNLmodel(photo=5, model="log",
       R1=R1, R2=Rb, Rb=Rb, I=D65, C=C,
       noise=TRUE, e = c(0.13, 0.06, 0.12, 0.07, 0.08))

#3
## Photoreceptor sensitivity spectra
##with lambda max at 350nm, 450nm and 550nm:
C<-photor(lambda.max=c(350,450,550))

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

## Run model 
model<-RNLmodel(photo=3, model="log",
       R1=R1, R2=R2, Rb=Rb, I=D65, C=C,
       noise=FALSE, v = c(NA, 0.06, NA),
       n = c(1,2,2))
}