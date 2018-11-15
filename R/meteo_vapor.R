#' Calculate meteorological variables related to water vapor

#' @filename meteo_vapor.R
#' @author Maik Renner, mrenner [at] bgc-jena.mpg.de
#' @export Magnus_Alduchov1996Improved
#' @export Magnus_Alduchov1996Improved_slope_sat

#' @version 0.10 2018-09-19 adding the Alduchov and Eskridge (1996) parameters for the Magnus form
#' @version 0.10 2018-09-19 adding an ice switch, be careful when vectorization is used for param ice
# weitere Sättigungsdampfdruckformeln unter http://cires.colorado.edu/~voemel/vp.html

Magnus_Alduchov1996Improved <- function(air_temp, c1 = 6.1094, c2 = 17.625, c3 = 243.04, ice = FALSE) {
  #' Empirical equation for saturation water vapor pressure

  #' based on the Magnus form
  #' coefficients taken from Alduchov and Eskridge (1996) Journal of Apllied Meteorology
  #'
  #' @param air_temp Temperature in °C.
  #' @param c1,c2,c3 Empirical coefficients with defaults set to values by
  #'        Alduchov and Eskridge (1996) Journal of Apllied Meteorology.
  #' @param ice Logical parameter which switches to the ice params
  #'  when an ice surface is there and T < 0 (default ice = FALSE).
  #' @return water vapor pressure at saturation over water surfaces or ice surface in hPa
  #' @details Alduchov and Eskridge (1996) Journal of Apllied Meteorology eqs 21 for water and 23 for ice
  #' @details also provide enhancement factors for moist air (not included here)
  #' @details see discussion by Koutsoyannis 2012 on the derivation
  #' @seealso [Magnus_Alduchov1996Improved_slope_sat]
  es = c1 * exp(c2 * air_temp / (c3 + air_temp))
## for ice surface
  if (ice == TRUE ) {
    c1 = 6.1121
    c2 = 22.587
    c3 = 273.86
    es[!is.na(air_temp) & air_temp < 0] = c1 * exp(c2 * air_temp[!is.na(air_temp) & air_temp < 0] / (c3 + air_temp[!is.na(air_temp) & air_temp < 0]))
  }
  return(es)
}

Magnus_Alduchov1996Improved_slope_sat <- function(air_temp,  c2 = 17.625, c3 = 243.04, ice = FALSE) {
  #' calc the derivative of the saturation vapor pressure curve of MAGNUS

  #' @param air_temp Temperature in °C.
  #' @param c2,c3 Empirical coefficients with defaults set to values by
  #'        Alduchov and Eskridge (1996) Journal of Apllied Meteorology.
  #' @param ice Logical parameter which switches to the ice params
  #'  when an ice surface is there and T < 0 (default ice = FALSE).
  #' @return slope of the saturation water vapor pressure curve hPa K-1.
  #' @seealso [Magnus_Alduchov1996Improved]
 
  s_airtemp <- c2 * c3 * MAGNUS(air_temp) / (c3 + air_temp)^2
  if (ice == TRUE) {
    c1 = 6.1121
    c2 = 22.587
    c3 = 273.86
    s_airtemp[!is.na(air_temp) & air_temp < 0] <- c2 * c3 * MAGNUS(air_temp[!is.na(air_temp) & air_temp < 0]) / (c3 + air_temp[!is.na(air_temp) & air_temp < 0])^2
  }
  return(s_airtemp)
}


MAGNUS <- function(air_temp) {
#' Empirical equation for saturation water vapor pressure
#'
#' @param air_temp Temperature in °C
#' @return water vapor pressure at saturation over water surfaces in hPa
#' @details separates vapor pressure over water and ice
#' @details wikipedia reports different coefficients
# oft werden aber auch andere gebraucht, abstimmung notwendig
  c1 = 6.1078 # in hPa
  c2 = 17.08085
  c3 = 234.175
  es = c1 * exp(c2 * air_temp / (c3 + air_temp))
#  if (air_temp < 0 ) {  # over ice
  c2 = 17.84362
  c3 = 245.425
  es[!is.na(air_temp) & air_temp < 0] = c1 * exp(c2 * air_temp[!is.na(air_temp) & air_temp < 0] / (c3 + air_temp[!is.na(air_temp) & air_temp < 0]))
  return(es)
}

## Steigung der S?ttigungsdampfdruckkurve
#  Ableitung der Magnusformel nach der Temperatur

slope_sat <- function(air_temp) {
  #' calc the derivative of the saturation vapor pressure curve of MAGNUS
  #' @update 2018-08-17 check for NA when in assignment
#source("MAGNUS.R")
  c2 = 17.08085
  c3 = 234.175
  s_airtemp <- c2 * c3 * MAGNUS(air_temp) / (c3 + air_temp)^2
#  if (air_temp < 0 ) {  # over ice
  c2 = 17.84362
  c3 = 245.425
  # s_airtemp[air_temp < 0] <- c2 * c3 * MAGNUS(air_temp[air_temp < 0]) / (c3 + air_temp[air_temp < 0])^2
  s_airtemp[!is.na(air_temp) & air_temp < 0] <- c2 * c3 * MAGNUS(air_temp[!is.na(air_temp) & air_temp < 0]) / (c3 + air_temp[!is.na(air_temp) & air_temp < 0])^2

  return(s_airtemp)
}
