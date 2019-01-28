#' Calculate the phase lag between two variables
#' 
#' Calculate the phase lag between two time series in time units, 
#' extending the Camuffo-Bernardi (1980) regression model with a harmonic
#' analysis.
#' The functions accompany a scientific manuscript
#' published in the Journal Hydrology and Earth System Sciences
#' as Renner et al., 2019 "Using phase lags to evaluate model biases
#'  in simulating the diurnal cycle of evapotranspiration:
#'   a case study in Luxembourg, Hydrol. Earth Syst. Sci., 23, 515-535,
#'    https://doi.org/10.5194/hess-23-515-2019, 2019.
#' It also contains some statistical utilities 
#' as well as simple meteorological functions.


#' @filename phaselag.R
#' @version 0.01 20180704 start collection of functions from paper code
#' @version 0.02 20180730 provide an artificial example with predefined phase lag to illustrate / test method
#' @version 0.03 20180808 new camuffo_phaselag_time() full function to compute the phase lag between variable Y and X + dX in minutes plus statistics in a wide format
#' @version 0.1.3 20181114 fixed camuffo_phaselag_time() with @import data.table 
#' @version 0.1.4 20181221 bugfix in camuffo_phaselag_time()  error occured because slope1 and slope2 where not found 
#' @author Maik Renner, mrenner [at] bgc-jena.mpg.de

#' @import data.table
#' @export phaselag_time
#' @export camuffo_phaselag_time

#' @references please cite the following paper when using this code:
#' 
#' Renner, M., Brenner, C., Mallick, K., Wizemann, H.-D., Conte, L., Trebs, I., Wei, J., Wulfmeyer, V., Schulz, K., and Kleidon, A.: Using phase lags to evaluate model biases in simulating the diurnal cycle of evapotranspiration: a case study in Luxembourg, Hydrol. Earth Syst. Sci., 23, 515-535, \url{https://doi.org/10.5194/hess-23-515-2019}, 2019.
#' 
#' @details
#' The main function to compute the phase lag in time units is [camuffo_phaselag_time()].
#' Other important utilities are [mlm.output.statlong()] which reorganizes lm() objects
#'  in a long table format with statistic|value.
#' Meteorological relationshipts based on a bulk type representations are 
#' [aerodynamic_conductance_withcanopy_Thom1975()], [LatentHeatFlux_PenmanMonteith()], 
#' [Magnus_Alduchov1996Improved()], [Magnus_Alduchov1996Improved_slope_sat()]
#' 
#' @keywords internal
"_PACKAGE"
#> [1] "_PACKAGE"

#R
require(data.table)

phaselag_time = function(slope1, slope2, nday, timeunitperday) {
#' Calculate the phase lag in time units
#' 
#' Helper function which is called by camuffo_phaselag_time()
#' 
#' @param slope1 Linear regression slope between two variables X and Y.
#' @param slope2 2nd order regression slope between the variable Y and the time derivative of X dX/dt.
#' @param nday Number of measurements per day.
#' @param timeunitperday Number of timesteps per day.
#' @return phase lag in time units between two variables.
#' @author Maik Renner, mrenner [at] bgc-jena.mpg.de, Luigi Conte contributed the harmonic analysis
#' @seealso camuffo_phaselag_time
#' @references please cite the following paper when using this code
#' Renner, M., Brenner, C., Mallick, K., Wizemann, H.-D., Conte, L., Trebs, I., Wei, J., Wulfmeyer, V., Schulz, K., and Kleidon, A.: Using phase lags to evaluate model biases in simulating the diurnal cycle of evapotranspiration: a case study in Luxembourg, Hydrol. Earth Syst. Sci., 23, 515-535, https://doi.org/10.5194/hess-23-515-2019, 2019.


 atan( (-slope2*2*pi/nday)/slope1) * (timeunitperday) / (2 *pi)
}

camuffo_phaselag_time = function(Y,X,dX, ...) {
#' Main function to calc the phase lag between variabe Y and X
#'
#' Derive the phase lag of two variables.
#' This is based on a harmonic transformation of a multi-linear regression of
#'  the form
#' Y = aX + b dX/dt + c
#' proposed by Camuffo and Bernardi (1982).
#' Renner et al., 2019 adapted their version to calculate the time difference
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
#' Renner, M., Brenner, C., Mallick, K., Wizemann, H.-D., Conte, L., Trebs, I., Wei, J., Wulfmeyer, V., Schulz, K., and Kleidon, A.: Using phase lags to evaluate model biases in simulating the diurnal cycle of evapotranspiration: a case study in Luxembourg, Hydrol. Earth Syst. Sci., 23, 515-535, https://doi.org/10.5194/hess-23-515-2019, 2019.
#' @seealso [phaselag_time()] for the core function
#' @examples
#' # example with 1h timesteps in minutes as unit
#' library(data.table)
#' nday = 24
#' timeunitperday = 24 * 60#' time unit minutes
#' dt = data.table(time = seq(0.5,timeunitperday, by = timeunitperday/nday))
#' dt[ , Date := as.IDate("2007-07-07")]
#' dt[ , ptime := time * 2* pi / timeunitperday]
#' dt[ , ptimecos := cos(ptime  - pi) ]
#' dt[ , Rsd := ifelse(ptimecos<0,0,800 * ptimecos)]
#' dt
#' dt[ , LEnolag :=  Rsd/2 + 30 * rnorm(1) , by = time ]
#' # create a series with a lag of one hour
#' (lagLE = 60 *  2 * pi / timeunitperday)
#' dt[ , ptimecos_lag1h := cos(ptime  - pi - lagLE)  ]
#' dt[ , LElag1h :=  ifelse(ptimecos_lag1h<0,0, 400 * ptimecos_lag1h ) + 30 * rnorm(1)  , by = time ]
#' # estimate the first order difference from the reference series (Rsd)
#' dt[ ,dRsd := c(NA, diff(Rsd))]
#' # now call the function camuffo_phaselag_time() 
#' dt[ , camuffo_phaselag_time(Y = LElag1h, X = Rsd, dX = dRsd, dt = .SD, nday = nday, timeunitperday = timeunitperday)]
#' dt[ , camuffo_phaselag_time(Y = LEnolag, X = Rsd, dX = rep(NA,24), dt = .SD, nday = nday, timeunitperday = timeunitperday)]
#' # plot stuff   
#' plot(Rsd ~ time, data = dt)
#' lines(LEnolag ~ time, data = dt, col = 3)
#' lines(LElag1h ~ time, data = dt, col = 4)
#
#' plot(LEnolag ~ Rsd, data = dt, col = 3, type = "b")
#' lines(LElag1h ~ Rsd, data = dt, col = 4, type = "b", pch = 0)
  
  # reg = mlm.output.statlong.call("Y ~ X + dX", dt)
  if (inherits(try(ans<-lm(Y ~ X + dX) ,silent = TRUE),"try-error")) {
    # reg = data.table(statistic = NA_character_, value = NA_real_)
    regwide = data.table()
   } else if (any(is.na(coef(ans)))) {
     # @version 2018-12-21 bugfix  error occured because slope1 and slope2 where not found 
     regwide = data.table()
    # data.table(statistic = NA_character_, value = NA_real_)
  } 
  else {
    reg = mlm.output.statlong(ans)
    regwide = dcast.data.table(reg, ... ~ statistic)
    regwide[ , '.' := NULL]
    if (reg[ , "slope1" %in% statistic] ) { 
        regwide[ , phaselagtime := phaselag_time(slope1, slope2, ... )]
    }
    return(data.table(regwide))
  }
}


#' @examples
#' ## Artificial example with solar radiation and timesteps in minutes as unit
#' 
#' # library(devtools)
#' # install_github("laubblatt/phaselag")
#' library(phaselag)
#' 
#' # Set number of time steps per period (day)
#' nday = 24
#' #' set time unit minutes
#' timeunitperday = 24 * 60
#' dt = data.table(time = seq(0.5,timeunitperday, by = timeunitperday/nday))
#' dt[ , Date := as.IDate("2007-07-07")]
#' dt[ , ptime := time * 2* pi / timeunitperday]
#' dt[ , ptimecos := cos(ptime  - pi) ]
#' dt[ , Rsd := ifelse(ptimecos<0,0,800 * ptimecos)]
#' dt
#' # create a second time series without a time lag
#' dt[ , LEnolag :=  Rsd/2 + 30 * rnorm(1) , by = time ]
#' # create another series which has a lag of one hour
#' (lagLE = 60 *  2 * pi / timeunitperday)
#' dt[ , ptimecos_lag1h := cos(ptime  - pi - lagLE)  ]
#' dt[ , LElag1h :=  ifelse(ptimecos_lag1h<0,0, 400 * ptimecos_lag1h ) + 30 * rnorm(1)  , by = time ]
#' 
#' plot(Rsd ~ time, data = dt)
#' lines(LEnolag ~ time, data = dt, col = 3)
#' lines(LElag1h ~ time, data = dt, col = 4)
#' 
#' ## Hysteresis loops appear when plotted against Reference
#' plot(LEnolag ~ Rsd, data = dt, col = 3, type = "b")
#' lines(LElag1h ~ Rsd, data = dt, col = 4, type = "b", pch = 0)
#' 
#' #### do regression and estimate time lags
#' # first estimate the time derivative of the reference series
#' dt[ ,dRsd := c(NA, diff(Rsd))]
#' # then call the estimation of the phase lags
#' dt[ , camuffo_phaselag_time(Y = LElag1h, X = Rsd, dX = dRsd, nday = nday, timeunitperday = timeunitperday)]
#' dt[ , camuffo_phaselag_time(Y = LEnolag, X = Rsd, dX = dRsd, nday = nday, timeunitperday = timeunitperday)]
#' ## compare  phaselagtime with defined phase lag
#' ## the significance of the lag can be assessed by slope2_pvalue
#' 
#' ## check error handling with missing input data 
#' dt[ , camuffo_phaselag_time(Y = LEnolag, X = Rsd, dX = rep(NA,24), nday = nday, timeunitperday = timeunitperday)]
#' #
#' 
#' ## Example which does regression and timelag estimation separately 
#' (camuffo_nolag =  dt[ , as.list(coef(lm(LEnolag ~ Rsd + dRsd))), by = Date])
#' (camuffo_lag1h =  dt[ , as.list(coef(lm(LElag1h ~ Rsd + dRsd))), by = Date])
#' ### compute the phase lag in time units (here time)
#' phaselag_time(slope1 = camuffo_nolag[ ,Rsd], slope2 = camuffo_nolag[ , dRsd], nday = nday, timeunitperday = timeunitperday)
#' phaselag_time(slope1 = camuffo_lag1h[ ,Rsd], slope2 = camuffo_lag1h[ , dRsd], nday = nday, timeunitperday = timeunitperday)


