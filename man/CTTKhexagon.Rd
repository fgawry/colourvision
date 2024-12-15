\name{CTTKhexagon}
\alias{CTTKhexagon}
\title{
Chittka (1992) colour hexagon}
\description{
Plots Chittka (1992) colour hexagon for trichromatic animals and a line plot for dichromatic animals.
}
\usage{
CTTKhexagon(x, y, photo=3,
            vnames=c(expression(E[1]),expression(E[2]),expression(E[3])),
            pch=16, bty="n", yaxt="n",xaxt="n", col="black",
            xlim="auto", ylim="auto", asp=1, ann=FALSE,
            axes=FALSE, vectors=FALSE, ...)
}
\arguments{
  \item{x}{x-coordinate}
  \item{y}{y-coordinate. \code{y=0} when \code{photo=2}}
  \item{photo}{Number of photoreceptor types. \code{photo=3} for a trichromatic animal, and \code{photo=2} for a dichromatic animal.}
  \item{vnames}{Vector names.}
  \item{pch}{see \code{\link{par}} function.}
  \item{bty}{see \code{\link{par}} function.}
  \item{yaxt}{see \code{\link{par}} function.}
  \item{xaxt}{see \code{\link{par}} function.}
  \item{col}{see \code{\link{par}} function.}
  \item{xlim}{see \code{\link{plot}} function. Default calculates xlim automatically.}
  \item{ylim}{see \code{\link{plot}} function. Default calculates ylim automatically.}
  \item{asp}{see \code{\link{plot}} function.}
  \item{ann}{see \code{\link{par}} function.}
  \item{axes}{see \code{\link{plot.default}} function.}
  \item{vectors}{Whether vectors representing direction of photoreceptor outputs should be plotted.}
  \item{...}{Other arguments passed to \code{\link{plot}} function.}
}
\references{
Chittka, L. 1992. The colour hexagon: a chromaticity diagram based on photoreceptor excitations as a generalized representation of colour opponency. J Comp Physiol A 170:533-543.
}
\author{
Felipe M. Gawryszewski \email{f.gawry@gmail.com}
}
\seealso{
\code{\link{CTTKmodel}}, \code{\link{CTTKhexagon3D}}
}
\examples{
##Honeybee photoreceptor sensitivity curves
data("bee")

##Grey background:
## with 7 percent reflectance from 300 to 700nm:
Rb <- data.frame(300:700, rep(7, length(300:700)))

## Read CIE D65 standard illuminant already converted to quantum flux:
data("D65")

##Reflectance data
## with a sigmoid spectrum and midpoint at 500nm and 550 nm
R1<-logistic(x=seq(300,700,1), x0=500, L=50, k=0.04)
R2<-logistic(x=seq(300,700,1), x0=550, L=50, k=0.04)
R<-cbind(R1, R2[,2])

## Run colour vision model:
model<-CTTKmodel(photo="tri", R=R, I=D65, Rb=Rb,
C=bee)

##plot data in the colour space
CTTKhexagon(x=model[,"X1"], y=model[,"X2"])
}