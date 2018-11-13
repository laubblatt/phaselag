#' collection of specific useful summary functions
#' @filename utils.aggregate.R

#' @author Maik Renner mrenner@bgc-jena.mpg.de
#' @version 1.00 2018-09-21 copied from "summary_fun.R"
#' @version 1.01 2018-09-21 improve documentation

#' @export meana
#' @export meann

#' @version 0.1 2017-09-22 added weighted unbiased variance, weighted.var()

meana <- function(x,na.percentage = 0.5,FUN=mean) {
  #' a function to aggregate data but set a limit to the data at least available
  #' 
  #' @param na.percentage maximal ratio of missing values 
  #' @author Maik Renner mrenner@bgc-jena.mpg.de
  if (sum(is.na(x))/length(x) <= na.percentage) {
    #     eval(parse(text = paste(FUN,"(x,na.rm = TRUE)",sep="")))
    eval(FUN(x,na.rm = TRUE))

 } else NA_real_
# tt = c(2,3,4,5,6,NA)
# meana(tt)
# meana(tt,1/6)
}


# http://stackoverflow.com/questions/12125364/why-does-median-trip-up-data-table-integer-versus-double
meann = function(x,nmin,...) {
#' an average mean which requires a minimum of data points nmin

#'  needed if missing data is missing :)
#' @param x vector of numeric values
#' @param nmin minimum of data points accepted to calculate the mean (integer)
#' @param ... further arguments to mean() remember especially to set na.rm=TRUE
#' @version 20141208 add non-missing data lower nmin
#' @author Maik Renner mrenner@bgc-jena.mpg.de
#' @examples
#' x = 1:10
#' meann(x,nmin = 11)
#' meann(NULL,nmin = 11)
#' x = c(1,2,NA,Inf,NaN)
#' meann(x,2,na.rm=TRUE)
#' x = c(1,2,NA,NA,NaN)
#' meann(x,2,na.rm=TRUE)
#' meann(x,3,na.rm=TRUE)
#' sum(!is.na(x))

  n = length(x)
  nnona = sum(!is.na(x))
 if (n < nmin | nnona < nmin) out = NA
 else out = mean(x,...)
 return(as.double(out))  ### data.table needs consitent output ... for data allocation
}

#' simple helper fun to calc an anomaly given the full data
anomin = function(x, ...) x - min(x,...)
anomean = function(x, ...) x - mean(x,...)



summaryfun_nmin = function(x,nmin,FUN = mean,...) {
#' a summary function applied to a vector which requires a minimum of data points nmin
#' 
#' needed if missing data is missing :)
#' add non-missing data lower nmin @20141208
#' @version 2015-09-16
#' @author Maik Renner mrenner@bgc-jena.mpg.de
  n = length(x)
  nnona = sum(!is.na(x))
 if (n < nmin | nnona < nmin) out = NA
 else out = eval(FUN(x,...))
 return(as.double(out))  ### data.table needs consitent output ... for data allocation
}
# @examples
# summaryfun_nmin(c(1,2,NA) , nmin = 2, FUN = mean, na.rm = TRUE)
# summaryfun_nmin(c(1,2,NA) , nmin = 2, FUN = max, na.rm = TRUE)
# summaryfun_nmin(c(2,NA) , nmin = 2, FUN = max, na.rm = TRUE)



rmse <- function(obs, pred) sqrt(mean((obs-pred)^2,na.rm=TRUE))

NSE = function(o,p) {
  #' computes the nash sutcliffe efficiency
  ok <- complete.cases(o, p)
  p = p[ok]
  o = o[ok]
  EF <- 1 - sum((p - o)^2)/sum((o - mean(o))^2)
  EF
}

meandifference = function(o,p) {
  ok <- complete.cases(o, p)
  p = p[ok]
  o = o[ok]
  mean(p-o)
}

SumSquaredErrors = function(o,p) {
  #' computes the sum of squared errors
  ok <- complete.cases(o, p)
  p = p[ok]
  o = o[ok]
  sum( (p-o)^2) / var(o)
}


## plot only within the quantile range given by a
plot95 = function(x,y,a=0.05,...) {
  xlim = quantile(x,c(a,1-a),na.rm=TRUE)
  ylim = quantile(y,c(a,1-a),na.rm=TRUE)
  plot(x,y,xlim = xlim, ylim = ylim, ...)
}

weighted.var <- function(x,wt) {
#' computes an unbiased weigthed variance using reliability weights
#' @version 0.1 copied and adapted from fun SDMtools::wt.var
        s = which(is.finite(x + wt)); wt = wt[s]; x = x[s] #remove NA info
	sumwt = sum(wt)
        xbar = sum(wt * x)/sumwt # calc weighted mean
        # WIKI Reliability weights If the weights are instead non-random (reliability weights), we can determine a correction factor to yield an unbiased estimator. Taking expectations we have,
#        https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Weighted_sample_variance
        wvar = sum(wt *(x-xbar)^2)* (sumwt /(sumwt^2 - sum(wt^2)))
	return(wvar) #return the variance
}
