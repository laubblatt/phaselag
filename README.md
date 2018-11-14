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
