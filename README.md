# R package phaselag
    
This repository provides an implemenation in the R programming language.  
    The package accompanies a scientific manuscript
    submitted to Hydrology and Earth System Sciences
    as Renner et al., 2018 "Estimating and understanding model bias
    in simulating the diurnal cycle of evapotranspiration:
    a case study in Luxembourg" available at https://www.hydrol-earth-syst-sci-discuss.net/hess-2018-310/ .
    Main functionality is to calculate the phase lag between two time series in time units, 
    extending the Camuffo-Bernardi (1980) regression model with a harmonic
    analysis.
    The package also contains some statistical utilities 
    as well as simple meteorological functions.
    
### To install the R package use:
```R
library(devtools)
install_github("laubblatt/phaselag")
 ```

### Small example with an artifical diurnal cycle 
First prepare the data. Here the reference is Rsd and there are two series 
generated which have no phase lag and one with a phase lag of one our. 
Both are contaminated with white noise. 
```R
library(phaselag)

# Set number of time steps per period (day)
nday = 24
#' set time unit minutes
 timeunitperday = 24 * 60
dt = data.table(time = seq(0.5,timeunitperday, by = timeunitperday/nday))
dt[ , Date := as.IDate("2007-07-07")]
dt[ , ptime := time * 2* pi / timeunitperday]
dt[ , ptimecos := cos(ptime  - pi) ]
dt[ , Rsd := ifelse(ptimecos<0,0,800 * ptimecos)]
dt
# create a second time series without a time lag
dt[ , LEnolag :=  Rsd/2 + 30 * rnorm(1) , by = time ]
# create another series which has a lag of one hour
(lagLE = 60 *  2 * pi / timeunitperday)
dt[ , ptimecos_lag1h := cos(ptime  - pi - lagLE)  ]
dt[ , LElag1h :=  ifelse(ptimecos_lag1h<0,0, 400 * ptimecos_lag1h ) + 30 * rnorm(1)  , by = time ]
 ```
Plot the generated series 
```R
plot(Rsd ~ time, data = dt)
lines(LEnolag ~ time, data = dt, col = 3)
lines(LElag1h ~ time, data = dt, col = 4)

## Hysteresis loops appear when plotted against Reference
plot(LEnolag ~ Rsd, data = dt, col = 3, type = "b")
lines(LElag1h ~ Rsd, data = dt, col = 4, type = "b", pch = 0)
 ```
Then perform regression and estimate time lags
```R
# first estimate the time derivative of the reference series
dt[ ,dRsd := c(NA, diff(Rsd))]
# then call the estimation of the phase lags
dt[ , camuffo_phaselag_time(Y = LElag1h, X = Rsd, dX = dRsd, nday = nday, timeunitperday = timeunitperday)]
dt[ , camuffo_phaselag_time(Y = LEnolag, X = Rsd, dX = dRsd, nday = nday, timeunitperday = timeunitperday)]
## compare  phaselagtime with defined phase lag
## the significance of the lag can be assessed by slope2_pvalue
 ```

Check error handling with missing input data 
```R
dt[ , camuffo_phaselag_time(Y = LEnolag, X = Rsd, dX = rep(NA,24), nday = nday, timeunitperday = timeunitperday)]
 ```

Example which does regression and timelag estimation separately 
```R
(camuffo_nolag =  dt[ , as.list(coef(lm(LEnolag ~ Rsd + dRsd))), by = Date])
(camuffo_lag1h =  dt[ , as.list(coef(lm(LElag1h ~ Rsd + dRsd))), by = Date])
### compute the phase lag in time units (here time)
phaselag_time(slope1 = camuffo_nolag[ ,Rsd], slope2 = camuffo_nolag[ , dRsd], nday = nday, timeunitperday = timeunitperday)
phaselag_time(slope1 = camuffo_lag1h[ ,Rsd], slope2 = camuffo_lag1h[ , dRsd], nday = nday, timeunitperday = timeunitperday)
 ```

