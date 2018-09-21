% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phaselag.R
\name{camuffo_phaselag_time}
\alias{camuffo_phaselag_time}
\title{Main function to calc the phase lag between variabe Y and X}
\usage{
camuffo_phaselag_time(Y, X, dX, ...)
}
\arguments{
\item{Y}{Time series of dependend variable.}

\item{X}{Time series of the reference variable, which should be close to a
harmonic.}

\item{dX}{Time derivative of X, eg. simply use c(NA,diff(X)).}

\item{...}{Further arguments to phaselag_time, set nday and timeunitperday.}
}
\value{
Returns a data.table with columns statistic and value
}
\description{
Derive the phase lag of two variables.
This is based on a harmonic transformation of a multi-linear regression of
 the form
Y = aX + b dX/dt + c
proposed by Camuffo and Bernardi (1982).
Renner et al., 2018 adapted their version to calculate the time difference
using the coefficients a and b determined by linear regression.
}
\details{
This function performs a multi-linear regression and then calls the function
phaselag_time() to estimate the time difference.
This function then also returns the regression statistics using the
function mlm.output.statlong()
}
\references{
please cite the following paper when using this code
Renner, M., Brenner, C., Mallick, K., Wizemann, H.-D., Conte, L., Trebs, I., Wei, J., Wulfmeyer, V., Schulz, K., and Kleidon, A.: Understanding model biases in the diurnal cycle of evapotranspiration: a case study in Luxembourg, Hydrol. Earth Syst. Sci. Discuss., https://doi.org/10.5194/hess-2018-310, in review, 2018.
}
\seealso{
[phaselag_time()]
}