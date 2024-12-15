\name{noise_e}
\alias{noise_e}
\title{
Receptor noise}
\description{
Receptor noise either provided by the user or based on noise-to-signal ratio of a single photoreceptor and the relative abundance of photoreceptor types in the retina. This function is used internally in Receptor Noise Limited models.
}
\usage{
noise_e(noise, e, v, n)
}
\arguments{
  \item{noise}{
Logical. Whether receptor noise is provided (\code{noise = TRUE}) or calculated from photoreceptor relative abundances (\code{noise = FALSE})
}
  \item{e}{
Receptor noise of each photoreceptor type. It is used when \code{noise = TRUE}
}
  \item{v}{
Noise-to-signal ratio of a single photoreceptor.
}
  \item{n}{
Relative abundance of the photoreceptor in the retina.
}
}
\references{
Vorobyev, M., and D. Osorio. 1998. Receptor noise as a determinant of colour thresholds. Proceedings of the Royal Society B 265:351-358.
}
\author{
Felipe M. Gawryszewski \email{f.gawry@gmail.com}
}
\seealso{
\code{\link{RNLmodel}}, \code{\link{RNLthres}}, \code{\link{GENmodel}}, \code{\link{colour_space}}
}