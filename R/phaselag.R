#' compute phase lag using Camuffo-Bernardi 1982 regression
#' adapted by Renner et al., 2018

#' @filename phaselag.R
#' @version 0.01 20180704 start collection of functions from paper code
#' @version 0.02 20180730 provide an artificial example with predefined phase lag to illustrate / test method
#' @version 0.03 20180808 new camuffo_phaselag_time() full function to compute the phase lag between variable Y and X + dX in minutes plus statistics in a wide format
#' @author mrenner@bgc-jena.mpg.de

#' @import data.table
#' @export

#' @reference please cite the following paper when using this code
 # Renner, M., Brenner, C., Mallick, K., Wizemann, H.-D., Conte, L., Trebs, I., Wei, J., Wulfmeyer, V., Schulz, K., and Kleidon, A.: Understanding model biases in the diurnal cycle of evapotranspiration: a case study in Luxembourg, Hydrol. Earth Syst. Sci. Discuss., https://doi.org/10.5194/hess-2018-310, in review, 2018.

#R
require(data.table)

phaselag_time = function(slope1, slope2, nday, timeunitperday) {
#' Calculate the phase lag in time units

#' @param slope1 Linear regression slope between two variables X and Y.
#' @param slope2 2nd order regression slope between the variable Y and the time derivative of X dX/dt.
#' @param nday Number of measurements per day.
#' @param timeunitperday Number of timesteps per day.
#' @return phase lag in time units between two variables.
#' @reference please cite the following paper when using this code
 # Renner, M., Brenner, C., Mallick, K., Wizemann, H.-D., Conte, L., Trebs, I., Wei, J., Wulfmeyer, V., Schulz, K., and Kleidon, A.: Understanding model biases in the diurnal cycle of evapotranspiration: a case study in Luxembourg, Hydrol. Earth Syst. Sci. Discuss., https://doi.org/10.5194/hess-2018-310, in review, 2018.

 atan( (-slope2*2*pi/nday)/slope1) * (timeunitperday) / (2 *pi)
}

# # dtseb_camufforeg_dRsd_widestat[ , phaselagtime :=
# atan( (-slope2*2*pi/48)  /slope1) * (60 * 24) / (2 *pi)]

### artificial example


# pr = "~/Dissertation/maiphd/R/"
## this requires to load data.table.regression.fun.R
# source(paste(pr,"data.table.regression.fun.R",sep=""))

camuffo_phaselag_time = function(Y,X,dX, ...) {
#' Main function to calc the phase lag between variabe Y and X
#'
#' Derive the phase lag of two variables.
#' This is based on a harmonic transformation of a multi-linear regression of
#'  the form
#' Y = aX + b dX/dt + c
#' proposed by Camuffo and Bernardi (1982).
#' Renner et al., 2018 adapted their version to calculate the time difference
#' using the coefficients a and b determined by linear regression.
#'
#' This function performs a multi-linear regression and then calls the function
#' phaselag_time() to estimate the time difference.
#' This function then also returns the regression statistics using the
#' function mlm.output.statlong()
#' @param Y Time series of dependend variable.
#' @param X Time series of the reference variable, which should be close to a
#' harmonic.
#' @param dX Time derivative of X, eg. simply use c(NA,diff(X)).
#' @param ... Further arguments to phaselag_time, set nday and timeunitperday.
#' @return Returns a data.table with columns statistic and value
#' @references please cite the following paper when using this code
#' Renner, M., Brenner, C., Mallick, K., Wizemann, H.-D., Conte, L., Trebs, I., Wei, J., Wulfmeyer, V., Schulz, K., and Kleidon, A.: Understanding model biases in the diurnal cycle of evapotranspiration: a case study in Luxembourg, Hydrol. Earth Syst. Sci. Discuss., https://doi.org/10.5194/hess-2018-310, in review, 2018.
#' @seealso [phaselag_time()]

  # reg = mlm.output.statlong.call("Y ~ X + dX", dt)
  if (inherits(try(ans<-lm(Y ~ X + dX) ,silent = TRUE),"try-error")) {
    reg = data.table(statistic = NA_character_, value = NA_real_)
    regwide = data.table()
  } else {
    reg = mlm.output.statlong(ans)
    regwide = dcast.data.table(reg, ... ~ statistic)
    regwide[ , '.' := NULL]
    # print(regwide)
    regwide[ , phaselagtime := phaselag_time(slope1, slope2, ... )]
    return(data.table(regwide))
  }
}


#' @examples
## example with 1h timesteps in minutes as unit
#' library(data.table)
#' nday = 24
#' timeunitperday = 24 * 60#' time unit minutes
#' dt = data.table(time = seq(0.5,timeunitperday, by = timeunitperday/nday))
#' dt[ , Date := as.IDate("2007-07-07")]
#' dt[ , ptime := time * 2* pi / timeunitperday]
#' dt[ , ptimecos := cos(ptime  - pi) ]
#' dt[ , Rsd := ifelse(ptimecos<0,0,800 * ptimecos)]
#' dt
#' plot(Rsd ~ time, data = dt)
#
#' dt[ , LEnolag :=  Rsd/2 + 30 * rnorm(1) , by = time ]
#' # impose a lag of one hour
#' (lagLE = 60 *  2 * pi / timeunitperday)
#' dt[ , ptimecos_lag1h := cos(ptime  - pi - lagLE)  ]
#' dt[ , LElag1h :=  ifelse(ptimecos_lag1h<0,0, 400 * ptimecos_lag1h ) + 30 * rnorm(1)  , by = time ]
#
#' lines(LEnolag ~ time, data = dt, col = 3)
#' lines(LElag1h ~ time, data = dt, col = 4)
#
#' plot(LEnolag ~ Rsd, data = dt, col = 3, type = "b")
#' lines(LElag1h ~ Rsd, data = dt, col = 4, type = "b", pch = 0)
#
#' ### do regression
#' dt[ ,dRsd := c(NA, diff(Rsd))]
#
#
#' (camuffo_nolag =  dt[ , as.list(coef(lm(LEnolag ~ Rsd + dRsd))), by = Date])
#
#' (camuffo_lag1h =  dt[ , as.list(coef(lm(LElag1h ~ Rsd + dRsd))), by = Date])
#
#' ## compute the phase lag in time units (here time)
#' phaselag_time(slope1 = camuffo_nolag[ ,Rsd], slope2 = camuffo_nolag[ , dRsd], nday = nday, timeunitperday = timeunitperday)
#
#' phaselag_time(slope1 = camuffo_lag1h[ ,Rsd], slope2 = camuffo_lag1h[ , dRsd], nday = nday, timeunitperday = timeunitperday)
#
#' mlm.output.statlong(lm(dt$LEnolag ~ dt$Rsd + dt$dRsd))
#
#' dt
#' colnames(dt)
#' dt[ , camuffo_phaselag_time(Y = LElag1h, X = Rsd, dX = dRsd, dt = .SD, nday = nday, timeunitperday = timeunitperday)]
#' dt[ , camuffo_phaselag_time(Y = LEnolag, X = Rsd, dX = rep(NA,24), dt = .SD, nday = nday, timeunitperday = timeunitperday)]
#