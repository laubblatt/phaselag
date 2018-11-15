#' Calculate different bulk conductance estimates

#' @filename bulk_conductance.R
#' @version 0.1.3 add functions to export to namespace 
#' @author Maik Renner, mrenner [at] bgc-jena.mpg.de

#' @export aerodynamic_conductance_withcanopy_Thom1975 aerodynamic_conductance_BM2013 aerodynamic_conductance_ustaru aerodynamic_conductance_Hinvert LatentHeatFlux_PenmanMonteith LatentHeatFlux_PenmanMonteith_gsinvert
#'

            
aerodynamic_conductance_BM2013 = function(ustar, u, karman = 0.4, Sc = 1, Pr = 1) {
  #' Calculate the aerodynamic conductance for heat from wind speed and friction velocity
  #' 
  #' Equation proposed by Baldocchi and Ma (2013) BLM
  #' 
  #' By using the measured friction velocity (u ∗ ) and wind speed (u) at the EC towers and using the equation of Baldocchi and Ma (2013) (g A-BM13 ) in which gA was expressed as the sum of turbulent conductance and canopy (quasi-laminar) boundary-layer conductance as
  #' result is very similar to ustar^2/u
  #' See also Brutsaert 1984 Evaporation into the atmosphere
  #' refferred in Verma 1989
  #' @param ustar friction velocity in m/s
  #' @param u horizontal wind speed in m/s
  #' @param karman Von Karmans' constant
  #' @param Sc Schmidt number 
  #' @param Pr Prandtl number 
  #' @return aerodynamic conductance in m/s
  #' @author Maik Renner, mrenner [at] bgc-jena.mpg.de
  ga_BM13 =  ( (u / ustar^2) + (2/ karman * ustar^2) * (Sc/Pr)^(2/3) )^-1
  return(ga_BM13)
}

aerodynamic_conductance_withcanopy_Thom1975 = function(ustar, u) {
  #' Empirical estimate of the aerodynamic conductance for heat from wind speed and friction velocity
  #' 
  #' Combined aerodynamic conductance after Thom 1975 taken from De Kauwe 2017 BG , eq. 3
  #' first term is the turbulent aerodynamic resistance the second term is the canopy boundary layer component
  #'
  #' @param ustar friction velocity in m/s
  #' @param u horizontal wind speed in m/s
  #' @return aerodynamic conductance in m/s
  #' @author Maik Renner, mrenner [at] bgc-jena.mpg.de
  
  
    ga =  ( (u / ustar^2) +  (6.2 * ustar^(-2/3)) )^-1
    return(ga)
}

aerodynamic_conductance_ustaru = function(ustar, u) { 
  #' Estimate of the aerodynamic conductance for momentum from wind speed and friction velocity
  #' @param ustar friction velocity in m/s
  #' @param u horizontal wind speed in m/s
  #' @return aerodynamic conductance in m/s

  ustar^2/u
}


aerodynamic_conductance_Hinvert = function(H, Ts, Ta, rho = 1.2, cp = 1004) {
  #' Estimate of the aerodynamic conductance for heat by inverting the bulk conductance formula for sensible heat
  #' 
  #' @param H sensible heat flux in W m-2
  #' @param Ts surface temperature in K, ideally this is the surface temperature of the heat source 
  #' @param Ta air temperature in K 
  #' @param rho density of air in kg m-3
  #' @param cp specific heat of air at constant pressure in J kg-1 K-1  
  #' @return aerodynamic conductance in m/s
  #' @author Maik Renner, mrenner [at] bgc-jena.mpg.de
  
  H / (rho * cp * (Ts - Ta))
}

LatentHeatFlux_PenmanMonteith = function(AE, s, esurf, eair, ga, gs, rho = 1.2, cp = 1004, psychro = 0.65) {
#' Calculate the Latent Heat flux using the Penman-Monteith formulation
#' 
#' @param AE Available Energy usually Net Radiation - Soil Heat flux in W m-2
#' @param s  slope of the saturation vapor pressure curve in hPa/K
#' @param esurf surface water vapour pressure hPa
#' @param eair vapour pressure in air hPa
#' @param ga aerodynamic conductance for water vapor in m/s 
#' @param gs surface / or stomatal conductance for water vapor in m/s 
#' @param rho Air density kg/m3
#' @param cp heat capacity of air J/(kg K)
#' @param psychro psychrometric constant hPa/K
#' @return Latent Heat flux in W/m2
#' @author Maik Renner, mrenner [at] bgc-jena.mpg.de
  
#' @seealso [aerodynamic_conductance_withcanopy_Thom1975()] ga needs a leaf boundary layer conductance (also known as excess resistance) 

(s * AE + rho * cp * ga  * (esurf - eair) ) / ( s + psychro* (1 + ga/gs) )

}

LatentHeatFlux_PenmanMonteith_gsinvert = function(AE, s, esurf, eair, ga, LE, rho = 1.2, cp = 1004, psychro = 0.65) {
  #' Calculate the surface conductance by inverting the Penman-Monteith formulation
  #' 
  #' @param AE Available Energy usually Net Radiation - Soil Heat flux in W m-2
  #' @param s  slope of the saturation vapor pressure curve in hPa/K
  #' @param esurf surface water vapour pressure hPa
  #' @param eair vapour pressure in air hPa
  #' @param ga aerodynamic conductance for water vapor in m/s 
  #' @param LE Latent Heat flux in W/m2
  #' @param rho Air density kg/m3
  #' @param cp heat capacity of air J/(kg K)
  #' @param psychro psychrometric constant hPa/K
  #' @return gs surface / or stomatal conductance for water vapor in m/s 
  #' @author Maik Renner, mrenner [at] bgc-jena.mpg.de
  #' @seealso [aerodynamic_conductance_withcanopy_Thom1975()] ga needs a leaf boundary layer conductance (also known as excess resistance) 
  gs = (LE*psychro*ga/((s*AE) - (LE*(s + psychro)) + (ga*rho*cp*(esurf - eair))))
}
