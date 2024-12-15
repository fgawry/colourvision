\name{plot.colourvision}
\alias{plot.colourvision}
\title{
Plot colour vision models into chromaticity diagrams
}
\description{
Plotting method for objects of class \code{colourvision}. Plotting method for animals with two or three photoreceptor types.
}
\usage{
\method{plot}{colourvision}(x, ...)
}

\arguments{
  \item{x}{
Object of class "\code{colourvision}".
}
  \item{\dots}{
Additional arguments passed to the plot function. See \code{\link{CTTKhexagon}} for Chittka (1992) model (\code{\link{CTTKmodel}}) plotting arguments; \code{\link{EMtriangle}} and \code{\link{EMline}} for Endler and Mielke (2005) model (\code{\link{EMmodel}}) plotting arguments; \code{\link{RNLplot}} for RNL plotting arguments; and \code{\link{GENplot}} for generic model (\code{\link{GENmodel}}) plotting arguments.
}
}
\references{
Chittka, L. 1992. The colour hexagon: a chromaticity diagram based on photoreceptor excitations as a generalized representation of colour opponency. J Comp Physiol A 170:533-543.

Endler, J. A., and P. Mielke. 2005. Comparing entire colour patterns as birds see them. Biol J Linn Soc 86:405-431.
}
\author{
Felipe M. Gawryszewski \email{f.gawry@gmail.com}
}
\seealso{
\code{\link{plot3d.colourvision}}, \code{\link{EMtriangle}}, \code{\link{CTTKhexagon}}, \code{\link{EMmodel}}, \code{\link{CTTKmodel}}, \code{\link{RNLmodel}}, \code{\link{RNLthres}}
}
\examples{
#trichromatic
##Photoreceptor sensitivity curves
C<-photor(lambda.max=c(350,450,550))

##Gray background
Rb <- data.frame(300:700, rep(7, length(300:700)))

## Read CIE D65 standard illuminant
data("D65")

##Reflectance data
R1<-logistic(x=seq(300,700,1), x0=500, L=50, k=0.04)
R1[,2]<-R1[,2]+10

##Run models
model<-EMmodel(photo=3,
       R=R1, I=D65, Rb=Rb, C=C)
plot(model)

model<-CTTKmodel(photo=3,
       R=R1, I=D65, Rb=Rb, C=C)
plot(model)

model<-RNLmodel(model="log", photo=3,
       R1=R1, I=D65, Rb=Rb, C=C, noise=TRUE, e=c(0.13, 0.06, 0.12))
plot(model)

#colour threshold
model<-RNLthres(photo=3, I=D65, Rb=Rb, C=C,
         noise=TRUE, e=c(0.13, 0.06, 0.12))
plot(model)

#dichromatic
##Photoreceptor sensitivity curves
C<-photor(lambda.max=c(400,550))

##Run models
model<-EMmodel(photo=2,
       R=R1, I=D65, Rb=Rb, C=C)
plot(model)

model<-EMmodel(photo=2, type="edge",
       R=R1, I=D65, Rb=Rb, C=C)
plot(model)

model<-CTTKmodel(photo=2,
       R=R1, I=D65, Rb=Rb, C=C)
plot(model)

model<-RNLmodel(model="log", photo=2,
       R1=R1, I=D65, Rb=Rb, C=C, noise=TRUE, e=c(0.13, 0.06))
plot(model)

#colour threshold
model<-RNLthres(photo=2, I=D65, Rb=Rb, C=C,
         noise=TRUE, e=c(0.13, 0.06))
plot(model)
}
