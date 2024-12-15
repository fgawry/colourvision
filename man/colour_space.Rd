\name{colour_space}
\alias{colour_space}
\title{
N-dimensional colour spaces
}
\description{
Generates a colour space based on any number of photoreceptor types and finds a colour locus for a given photoreceptor output.}
\usage{
colour_space(n, type="length", length=NA, edge=NA,
          q=rep(1,n), recep.noise=FALSE, e=NA)
}
\arguments{
  \item{n}{Number of photoreceptor types. Function accepts any number of photoreceptor types >=2. For instance, trichromatic: \code{n=3}; tetrachromatic: \code{n=4}, etc.
}
  \item{type}{Whether the colour space should be built with a fixed vector length (\code{type="length"}), or a fixed distance between vertices (\code{type="edge"}; e.g. Pike 2012 and Renoult et al. 2015).
}
  \item{length}{Vector length used to construct the colour space. Used when \code{type="length"}. Typically \code{length=1}.}
  \item{edge}{Edge length used to construct the colour space. Used when \code{type="edge"}.}
  \item{q}{Photoreceptor output values}
  \item{recep.noise}{Whether receptor noise should be used to calculate colour locus coordinates.}
  \item{e}{Vector representing photoreceptor noises. Used when \code{recep.noise=TRUE}}
}
\details{This function is used internally in colour vision models.
}
\value{
A list with the following dimensions:
\item{coordinates}{Colour locus of photoreceptor output values \code{q} in the colour space}
\item{vector_matrix}{Matrix of column vectors, each representing one photoreceptor type}
}
\references{
Pike, T.W. 2012. Generalised chromaticity diagrams for animals with n-chromatic colour vision. Journal of Insect Behavior 255: 277-286.

Renoult, J. P., A. Kelber, and H. M. Schaefer. 2015. Colour spaces in ecology and evolutionary biology. Biol Rev Camb Philos Soc, doi: 10.1111/brv.12230.
}
\author{
Felipe M. Gawryszewski \email{f.gawry@gmail.com}
}
\seealso{
\code{\link{Q}}, \code{\link{Qr}}, \code{\link{CTTKmodel}}, \code{\link{EMmodel}}, \code{\link{RNLmodel}}, \code{\link{GENmodel}}
}

\examples{
#A trichromatic colour space based on Endler and Mielke (2005) 
tri<-colour_space(n=3, length=0.75, q=c(0.5,0.2,0.3))

#showing:
#(1) Limits of the colour space (triangle)
plot(0, ylim=c(-1,1), xlim=c(-1,1), asp=1, ylab="X2", xlab="X1", type="n")
polygon(x=tri$vector_matrix[1,], y=tri$vector_matrix[2,], lty=2)

#(2) Vectors (length=0.75) used to build the colour space (arrows)
arrows(x0=0,y0=0, x1=tri$vector_matrix[1,1], y1=tri$vector_matrix[2,1], col="red")
arrows(x0=0,y0=0, x1=tri$vector_matrix[1,2], y1=tri$vector_matrix[2,2], col="red")
arrows(x0=0,y0=0, x1=tri$vector_matrix[1,3], y1=tri$vector_matrix[2,3], col="red")

#(3) Colour loci of given photoreceptor outputs
points(x=tri$coordinates[[1]], y=tri$coordinates[[2]], pch=21, col="blue", bg="blue")
}