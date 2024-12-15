\name{GENmodel}
\alias{GENmodel}
\title{
N-dimensional generic colour vision model
}
\description{
A flexible function to build colour vision models based on any number of photoreceptor types (Gawryszewski 2018).
}
\usage{
GENmodel(photo=ncol(C)-1, type="length", length=NA, edge=NA,
          R, I, Rb=NA, C, vonKries = TRUE, func, unity=FALSE,
          recep.noise=FALSE, noise.given=TRUE, e=NA, v=NA, n=NA,
          interpolate=TRUE, nm=seq(300,700,1))
}
\arguments{
  \item{photo}{
Number of photoreceptor types. Model accepts any number of photoreceptor types (\code{>=2}). For instance, dichromatic: \code{photo=2}; trichromatic: \code{photo=3}; tetrachromatic: \code{photo=4}, etc. Default gets number of photoreceptor types from \code{C} argument.
}
  \item{type}{Whether the colour space should be built with a fixed vector length (\code{type="length"}), or a fixed distance between vertices (\code{type="edge"}; e.g. Pike 2012 and Renoult et al. 2015).
}
  \item{length}{Vector length used to construct the colour space. Used when \code{type="length"}. Typically \code{length=1}.}
  
  \item{edge}{Edge length used to construct the colour space. Used when \code{type="edge"}.}
  
  \item{R}{
Reflectance of observed objects. A data frame with first column corresponding to wavelength values and following columns with reflectance values. \code{R} must be in the same scale as \code{Rb}.
}
  \item{I}{
Irradiance spectrum. A data frame with two columns only: first column corresponding to wavelength values and second column with irradiance values. Irradiance values must be in quantum flux units.
}
  \item{Rb}{
Background reflectance. A data frame with two columns only: first column corresponding to wavelength values and second column with reflectance values. \code{Rb} must be in the same scale as in \code{R}.
}
  \item{C}{
Photoreceptor sensitivity curves, from lowest to longest lambda-max. A data frame: first column corresponding to wavelength values and following columns with photoreceptor sensitivity values (see function \code{\link{photor}}). 
}
  \item{vonKries}{Whether the von Kries transformation should be applied to photoreceptor photon catches. Calculates the ratio between stimulus and background photoreceptor photon catches}

  \item{func}{Transformation of photoreceptors inputs. For instance, \code{func=log} for a ln-transformation.}
  
  \item{unity}{Whether sum of photoreceptor output should be normalized to 1, as in \code{\link{EMmodel}}}
  
  \item{recep.noise}{Logical. Whether chromaticity distances should be calculated based on receptor noises.}

  \item{noise.given}{Logical. Whether receptor noise is provided (\code{noise.given = TRUE}) or calculated from photoreceptor relative abundances (\code{noise.given = FALSE}). Used when \code{recep.noise=TRUE}.}

  \item{e}{
Receptor noise of each photoreceptor type. Used when \code{recep.noise = TRUE} and \code{noise.given = TRUE}
}
  \item{n}{
Relative number of each photoreceptor type in the retina. Used to calculate \code{e} when \code{recep.noise = TRUE} and \code{noise.given = TRUE}.
}
  \item{v}{
Noise-to-signal ratio of a single photoreceptor. Used to calculate \code{e} when \code{recep.noise = TRUE} and \code{noise.given = TRUE}.
}
  \item{interpolate}{
Whether data files should be interpolated before further calculations. See \code{\link{approx}}.}

  \item{nm}{
A sequence of numeric values specifying where interpolation is to take place. See \code{\link{approx}}.}
}

\value{
\item{ei}{Noise of photoreceptor channels. Calculated only when \code{recep.noise=TRUE}}
\item{Qri}{Photoreceptor photon catch values from \code{R}. Relative Photoreceptor photon catches when \code{vonKries=TRUE}. See functions \code{\link{Q}} and \code{\link{Qr}}.}
\item{Ei}{Photoreceptor outputs.}
\item{Xi}{Stimulus colour locus coordinates in the colour space}
\item{deltaS}{Euclidean distance to the origin of the colour space. It represents the conspicuousness of the stimulus (\code{R}) in relation to the background (\code{Rb})}
}
\references{
Gawryszewski, F.M. 2018. Colour vision models: Some simulations, a general n-dimensional model, and the colourvision R package. Ecology and Evolution, 10.1002/ece3.4288.

Pike, T.W. 2012. Generalised chromaticity diagrams for animals with n-chromatic colour vision. Journal of Insect Behavior 255: 277-286.

Renoult, J. P., A. Kelber, and H. M. Schaefer. 2015. Colour spaces in ecology and evolutionary biology. Biol Rev Camb Philos Soc, doi: 10.1111/brv.12230.
}
\author{
Felipe M. Gawryszewski \email{f.gawry@gmail.com}
}
\seealso{
\code{\link{Q}}, \code{\link{Qr}}, \code{\link{CTTKmodel}}, \code{\link{EMmodel}}, \code{\link{RNLmodel}}, \code{\link{colour_space}}
}

\examples{
#A trichromatic colour vision model based on Endler and Mielke (2005)

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
model<-GENmodel(length=0.75, R=R, I=D65, Rb=Rb, C=C, 
                func=log, unity=TRUE)

plot(model)
}