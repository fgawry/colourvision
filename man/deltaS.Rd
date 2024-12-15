\name{deltaS}
\alias{deltaS}
\title{
Chromaticity distances
}
\description{
Calculates a matrix with all possible pairwise comparison between stimulus reflectance spectra based on a given colour vision model output.}
\usage{
deltaS(model)
}
\arguments{
  \item{model}{Output of a colour vision model}.
}

\value{A matrix with pairwise chromaticity distances.}

\author{
Felipe M. Gawryszewski \email{f.gawry@gmail.com}
}
\seealso{
\code{\link{CTTKmodel}}, \code{\link{RNLmodel}}, \code{\link{EMmodel}},\code{\link{GENmodel}}
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

#Chromaticity distances between R1, R2 and R3
deltaS(model)
}

