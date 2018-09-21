#' Calculate different bulk conductance estimates

#' @filename bulk_conductance.R

##### aerodynamic and surface conductance
aerodynamic_conductance_BM2013 = function(ustar, u, karman = 0.4, Sc = 1, Pr = 1) {
    #' By using the measured friction velocity (u ∗ ) and wind speed (u) at the EC towers and using the equation of Baldocchi and Ma (2013) (g A-BM13 ) in which gA was expressed as the sum of turbulent conductance and canopy (quasi-laminar) boundary-layer conductance as
    #' result is very similar to ustar^2/u
    #' See also Brutsaert 1984 Evaporation into the atmosphere
    #' refferred in Verma 1989
    ga_BM13 =  ( (u / ustar^2) + (2/ karman * ustar^2) * (Sc/Pr)^(2/3) )^-1
    return(ga_BM13)
    }

aerodynamic_conductance_withcanopy_Thom1975 = function(ustar, u) {
    #' combined aerodynamic conductance after Thom 1975 taken from De Kauwe 2017 BG , eq. 3
    #' first term is the turbulent aerodynamic resistance the second term is the canopy boundary layer component
    #'
    ga =  ( (u / ustar^2) +  (6.2 * ustar^(-2/3)) )^-1
    return(ga)
}


aerodynamic_conductance_ustaru = function(ustar, u) ustar^2/u
aerodynamic_conductance_Hinvert = function(H, Ts, Ta, rho = 1.2, cp = 1004) H / (rho * cp * (Ts - Ta))

LatentHeatFlux_moisturetransfer = function(esurf, eair, ra, rs, rho = 1.2, cp = 1004, psychro = 0.65) {
#' @param rho Air density kg/m3
#' @param cp heat capacity of air J/(kg⋅K)
#' @param psychro psychrometric constant hPa/K
#' @param esurf surface water vapour pressure hPa; actually this should be e0 as we do not consider gs here !
#' @param ra aerodynamic resistance s/m
#' @param

rho * cp / psychro * (esurf - eair) / (ra + rs)
}

LatentHeatFlux_moisturetransfer_rsinvert = function(LE, esurf, eair, ra, rho = 1.2, cp = 1004, psychro = 0.65) {
#' @param rho Air density kg/m3
#' @param cp heat capacity of air J/(kg⋅K)
#' @param psychro psychrometric constant hPa/K
#' @param esurf surface water vapour pressure hPa
#' @param ra aerodynamic resistance s/m
#' @param LE LatentHeatFlux W/m2
rho * cp / psychro * (esurf - eair) / LE - ra
}


LatentHeatFlux_PenmanMonteith = function(AE, s, esurf, eair, ga, gs, rho = 1.2, cp = 1004, psychro = 0.65) {
#' @param rho Air density kg/m3
#' @param cp heat capacity of air J/(kg⋅K)
#' @param psychro psychrometric constant hPa/K
#' @param esurf surface water vapour pressure hPa
#' @param ra aerodynamic resistance s/m
#  @param ga needs a leaf boundary layer conductance (also excess resistance) as well use ga_thom
#' @param

(s * AE + rho * cp * ga  * (esurf - eair) ) / ( s + psychro* (1 + ga/gs) )
}

LatentHeatFlux_PenmanMonteith_gsinvert = function(AE, s, esurf, eair, ga, LE, rho = 1.2, cp = 1004, psychro = 0.65) {
#' @param rho Air density kg/m3
#' @param cp heat capacity of air J/(kg⋅K)
#' @param psychro psychrometric constant hPa/K
#' @param esurf surface water vapour pressure hPa
#' @param ra aerodynamic resistance s/m
#' @param
gs = (LE*psychro*ga/((s*AE) - (LE*(s + psychro)) + (ga*rho*cp*(esurf - eair))))
}
