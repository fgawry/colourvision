\name{radarplot}
\alias{radarplot}
\title{
Radar plot
}
\description{
Plots quantum catches or E-values (photoreceptor outputs) into a radar plot.}
\usage{
radarplot(model, item=c("Qr", "E"), item.labels=FALSE, item.lwd=1,
           border=NULL, radar.lwd=1, radar.col="grey",
           length="auto", xlim="auto", ylim="auto",
           xlab="", ylab="", asp=1, add = FALSE, ...)
}
\arguments{
  \item{model}{Output of a colour vision model.}
  \item{item}{Whether photoreceptor inputs (\code{Qr}) or outputs (\code{E}) should be plotted.}
  \item{item.labels}{Logical. Whether vector names should be plotted.}
  \item{item.lwd}{Width of lines connecting item values.}
  \item{border}{Colour of lines connecting item values. See \code{\link{polygon}}.}
  \item{radar.lwd}{Width of lines representing item vectors.}
  \item{radar.col}{Colour of lines representing item vectors.}
  \item{length}{Length of \code{item} vectors. Default calculates length automatically.}
  \item{xlim}{x-axis range. Default calculates range automatically. See \code{\link{plot}} function.}
  \item{ylim}{y-axis range. Default calculates range automatically. See \code{\link{plot}} function.}
  \item{xlab}{see \code{\link{par}} function.}
  \item{ylab}{see \code{\link{par}} function.}
  \item{asp}{see \code{\link{par}} function.}
  \item{add}{Logical. Whether data should be added to an existing plot.}
  \item{...}{Other arguments passed to \code{\link{plot}}.}   
}
\author{
Felipe M. Gawryszewski \email{f.gawry@gmail.com}
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
R1<-logistic(x=seq(300,700,1), x0=450, L=50, k=0.04)
R2<-logistic(x=seq(300,700,1), x0=500, L=50, k=0.04)
R3<-logistic(x=seq(300,700,1), x0=550, L=50, k=0.04)
R<-cbind(R1,R2[,2],R3[,2])

## Run model 
model<-CTTKmodel(photo=3, R=R, I=D65,
    Rb=Rb, C=C)

#Radarplot
radarplot(model, border=c("violet", "red", "blue"), item="E", item.labels=TRUE)
}