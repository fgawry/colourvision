\name{energytoflux}
\alias{energytoflux}
\title{
Irradiance from energy to quantum units.
}
\description{
Convert Irradiance datum from energy units (uW/cm2/nm) to quantum flux units (umol/m2/s)
}
\usage{
energytoflux(datum)
}
\arguments{
  \item{datum}{A data frame with first column corresponding to wavelength values and second column with irradiance values. Irradiance values must be in energy units (uW/cm2/nm).
}
}
\value{
A data frame with first column corresponding to wavelength values and second column with irradiance values in umol/m2/s.
}
\author{
Felipe M. Gawryszewski \email{f.gawry@gmail.com}
}