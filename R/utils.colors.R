#' small utilities functions to deal with colors in R

#' @filename utils.colors.R
#' @author Maik Renner, mrenner [at] bgc-jena.mpg.de
#' @version 1.00 2018-09-21 copied from "utils.color.R"
#' @version 1.01 2018-09-21 improve documentation adding example

#' @export makeTransparent


makeTransparent<-function(someColor, alpha=100)
{
#' calculates the color code with a desired transparancy of a given color

#' useful for drawing polygons

#' @param someColor any color code what R understands as character
#' @param alpha numeric 0 ... 100 setting the transparancy, 0 no color,
#'  100 full color
#' @return return a color code which can be used in plot(  , col = )
#' @author Maik Renner, mrenner [at] bgc-jena.mpg.de
#' @examples
#' makeTransparent("green")
#' makeTransparent("green", alpha = 0.1)
#' makeTransparent(makeTransparent("green", alpha = 0.1))
#' plot(0,0, pch = 16, col = "green", cex = 5)
#' points(0.1,0, pch = 15, col = makeTransparent("green", alpha = 50), cex = 5 )
#' set.seed(688)
#' xp = rnorm(7)
#' yp = rnorm(7)
#' polygon(xp, yp, border = 3, col = makeTransparent("green", alpha = 80) )

#' @seealso <https://stackoverflow.com/questions/8047668/transparent-equivalent-of-given-color>
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
    blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}
