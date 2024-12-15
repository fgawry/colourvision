\name{photor}
\alias{photor}
\title{
Photoreceptor sensitivity spectra.
}
\description{
Generates photoreceptor sensitivity spectra based on lambda-max values.
}
\usage{
photor(lambda.max, lambda = seq(300, 700, 1), beta.band = FALSE)
}
\arguments{
  \item{lambda.max}{
Vector with photoreceptor wavelength at maximum sensitivity values, in increasing order. 
}
  \item{lambda}{
Range and interval to calculate the sensitivity curves.
}
  \item{beta.band}{
Logical. Whether or not to include a beta-band in the sensitivity curve.}
}
\value{
A data frame with first column corresponding to wavelength values and following columns with photoreceptor sensitivity values
}
\references{
Govardovskii, V. I., N. Fyhrquist, T. Reuter, D. G. Kuzmin, and K. Donner. 2000. In search of the visual pigment template. Vis. Neurosci. 17:509-528.
}
\author{
Felipe M. Gawryszewski \email{f.gawry@gmail.com}
}
\seealso{
\code{\link{CTTKmodel}}, \code{\link{EMmodel}}, \code{\link{RNLmodel}}, \code{\link{RNLthres}}
}
\examples{
## Generates photoreceptor sensitivity
## values with lambda max at 350nm, 450nm and 550nm:
C<-photor(lambda.max=c(350,450,550))

plot(C[,2]~C[,1], type="l", col="violet")
lines(C[,3]~C[,1], type="l", col="blue")
lines(C[,4]~C[,1], type="l", col="green")
}