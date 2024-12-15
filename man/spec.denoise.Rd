\name{spec.denoise}
\alias{spec.denoise}
\title{
Smooth function for reflectance spectra.
}
\description{
Applies a \code{\link{smooth.spline}} for data frame containing spectrometric data.
}
\usage{
spec.denoise(specfiles, spar = 0.7, ...)
}
\arguments{
  \item{specfiles}{
A data frame with first column representing wavelength values and following columns with reflectance data.
}
  \item{spar}{see \code{\link{smooth.spline}}. Smoothing parameter, typically (but not necessarily) in (0,1].
}
  \item{...}{Other arguments passed to function \code{\link{smooth.spline}}.}

}
\value{
A data frame with first column representing wavelength values and following columns with reflectance data.}
\author{
Felipe M. Gawryszewski \email{f.gawry@gmail.com}
}