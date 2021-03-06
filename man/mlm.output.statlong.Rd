% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.table.regression.utils.R
\name{mlm.output.statlong}
\alias{mlm.output.statlong}
\title{Reorganize the output of a regression fit into a statistic | value format}
\usage{
mlm.output.statlong(fit)
}
\arguments{
\item{fit}{object returned from a linear model function lm() and quantreg::rq() are supported}
}
\value{
a data.table with columns statistic and value
}
\description{
simple function which takes a lm model output to create a named vector intended for use with data.table and group by
}
\examples{
# output of regression model is melted to two columns
data(mtcars)
fitlm = lm(formula = "mpg ~ drat + gear", data = mtcars)
mlm.output.statlong(fitlm)
# example for group by operation with data.table power
DT = as.data.table(mtcars)
DT[, mlm.output.statlong(lm(formula = "mpg ~ drat + gear", .SD) ), by = list(am) ]
}
\author{
Maik Renner, mrenner [at] bgc-jena.mpg.de
}
