% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phaselag.R
\docType{package}
\name{phaselag-package}
\alias{phaselag}
\alias{phaselag-package}
\title{Calculate the phase lag between two variables}
\description{
Calculate the phase lag between two time series in time units, 
extending the Camuffo-Bernardi (1980) regression model with a harmonic
analysis.
The functions accompany a scientific manuscript
published in the Journal Hydrology and Earth System Sciences
as Renner et al., 2019 "Using phase lags to evaluate model biases
 in simulating the diurnal cycle of evapotranspiration:
  a case study in Luxembourg, Hydrol. Earth Syst. Sci., 23, 515-535,
   https://doi.org/10.5194/hess-23-515-2019, 2019.
It also contains some statistical utilities 
as well as simple meteorological functions.
}
\details{
The main function to compute the phase lag in time units is [camuffo_phaselag_time()].
Other important utilities are [mlm.output.statlong()] which reorganizes lm() objects
 in a long table format with statistic|value.
Meteorological relationshipts based on a bulk type representations are 
[aerodynamic_conductance_withcanopy_Thom1975()], [LatentHeatFlux_PenmanMonteith()], 
[Magnus_Alduchov1996Improved()], [Magnus_Alduchov1996Improved_slope_sat()]
}
\references{
please cite the following paper when using this code:

Renner, M., Brenner, C., Mallick, K., Wizemann, H.-D., Conte, L., Trebs, I., Wei, J., Wulfmeyer, V., Schulz, K., and Kleidon, A.: Using phase lags to evaluate model biases in simulating the diurnal cycle of evapotranspiration: a case study in Luxembourg, Hydrol. Earth Syst. Sci., 23, 515-535, \url{https://doi.org/10.5194/hess-23-515-2019}, 2019.
}
\author{
Maik Renner, mrenner [at] bgc-jena.mpg.de
}
\keyword{internal}
