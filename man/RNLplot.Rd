\name{RNLplot}
\alias{RNLplot}
\title{
Receptor noise limited model 2D and 1D plot}
\description{
Plots receptor noise limited model (RNL) for trichromatic and dichromatic animals.
}
\usage{
RNLplot(model, photo, item="R1",
        vectors=TRUE, vnames=TRUE, vsize="auto",
        xlab="x", ylab="y", xlim="auto", ylim="auto", asp=1, ...)
}
\arguments{
  \item{model}{Output of a colour vision model.}
  \item{photo}{Number of photoreceptor types.}
  \item{item}{Model output item to be plotted. Default plots stimulus data. See \code{\link{RNLmodel}}.}
  \item{vectors}{Whether vectors representing direction of photoreceptor outputs should be plotted.}
  \item{vnames}{Whether vector names should be plotted.}
  \item{vsize}{Length of vectors to be plotted. Default calculates length automatically.}
  \item{xlab}{x-axis range. Default calculates range automatically. See \code{\link{par}} function.}
  \item{ylab}{y-axis range. Default calculates range automatically. See \code{\link{par}} function.}
  \item{xlim}{see \code{\link{par}} function. Default calculates xlim automatically.}
  \item{ylim}{see \code{\link{par}} function. Default calculates ylim automatically.}
  \item{asp}{see \code{\link{plot}} function.}
  \item{...}{Other arguments passed to \code{\link{plot}} function.}
}
\author{
Felipe M. Gawryszewski \email{f.gawry@gmail.com}
}
\seealso{
\code{\link{CTTKhexagon}}, \code{\link{CTTKhexagon3D}}, \code{\link{EMtriangle}}, \code{\link{EMtetrahedron}}, \code{\link{RNLplot3d}}, \code{\link{plot.colourvision}}, \code{\link{plot3d.colourvision}}
}
\examples{
#dichromat
C<-photor(lambda.max=c(450,550))
Rb <- data.frame(300:700, rep(7, length(300:700)))
data("D65")
R1.1<-logistic(x=seq(300,700,1), x0=500, L=50, k=0.04)
R1.2<-logistic(x=seq(300,700,1), x0=400, L=50, k=0.04)
w<-R1.1[,1]
R1.1<-R1.1[,2]+10
R1.2<-R1.2[,2]+10
R1<-data.frame(w=w, R1.1=R1.1, R1.2=R1.2)
model<-RNLmodel(model="log",
       R1=R1, Rb=Rb, I=D65, C=C,
       noise=TRUE, e = c(0.13, 0.06))
plot(model)

#trichromat
C<-photor(lambda.max=c(350,450,550))
Rb <- data.frame(300:700, rep(7, length(300:700)))
data("D65")
R1.1<-logistic(x=seq(300,700,1), x0=500, L=50, k=0.04)
R1.2<-logistic(x=seq(300,700,1), x0=400, L=50, k=0.04)
w<-R1.1[,1]
R1.1<-R1.1[,2]+10
R1.2<-R1.2[,2]+10
R1<-data.frame(w=w, R1.1=R1.1, R1.2=R1.2)
model<-RNLmodel(model="log",
       R1=R1, Rb=Rb, I=D65, C=C,
       noise=TRUE, e = c(0.13, 0.06, 0.12))
plot(model)
}