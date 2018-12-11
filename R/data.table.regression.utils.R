#' functions which aid cross variable linear regressions for use with data.table
#' 
#' @filename data.table.regression.utils.R
#' @author Maik Renner, mrenner [at] bgc-jena.mpg.de

#' @version 1.00 2018-09-21 copies functions mlm.output.statlong() and mlm.output.statlong.call() from data.table.regression.fun.R
#' @version 1.01 2018-09-21 improve documentation

#' @import data.table
#' @export mlm.output.statlong
#' @export mlm.output.statlong.call

require(data.table)

mlm.output.statlong = function(fit) {
#' Reorganize the output of a regression fit into a statistic | value format 
#' 
#' simple function which takes a lm model output to create a named vector intended for use with data.table and group by
#' @param fit object returned from a linear model function lm() and quantreg::rq() are supported
#' @return a data.table with columns statistic and value
#' @author Maik Renner, mrenner [at] bgc-jena.mpg.de
#' @examples
#' # output of regression model is melted to two columns
#' data(mtcars)
#' fitlm = lm(formula = "mpg ~ drat + gear", data = mtcars)
#' mlm.output.statlong(fitlm)
#' # example for group by operation with data.table power
#' DT = as.data.table(mtcars)
#' DT[, mlm.output.statlong(lm(formula = "mpg ~ drat + gear", .SD) ), by = list(am) ]

#' @version 0.1 2016-03-11
#' @version 0.11 2016-04-15 add n column
#' @version 0.20 2016-05-05 long output format
#' @version 0.21 @20160713 added a error checking when coefficients matrix has less than 2 rows
#' @version 0.22 2017-03-16  adapted to different output of the summary()$coefficients object wher quantreg::rq had 3 column
#' @version 0.23 2018-02-02 rename output from variable to statistic
#' @version 0.24 2018-05-04 make the case for intercept == 0
#' @version 0.25 2018-08-06 error handling when summary(fit) throws an error
#' @version 0.25 2018-08-06 if fit is possible but summary(fit) throws an error (quantreg::qr ERROR singular matrix in 'backsolve') then a new if statement is triggered without error margins
#' @version 0.25 2018-08-06 also R2 is computed for all cases now
#' @version 0.25 2018-08-08 fit$y is only returned by rq() not by lm() so I used the more generic fit$model[,1] to access the response variable
#' @version 0.26 2018-11-   addding a call fun for quantreg 
#' @version 0.26 2018-11-29 removing a bug when regression forced through origin 
#' @version 0.26 2018-11-29 adding the case when coefficients are only a named num 
#' @version 0.26 2018-11-29 TODO check if regression still works with phaselag then add a pull request 
#' 
#' 

  #  print(mean(wght))
#  fit = lm(formula = frm, data = dt, weights = eval(wght))
 #fit = lm(frm, dt,weights = wght)
#  test = data.frame(Y = rnorm(90), X1 = rnorm(90), X2 = c(NA,NA,rnorm(88)) )
#  fit = lm(Y  ~ X1 + X2, data = test)

if (! is.null(fit) ) {
  if (inherits(try(sufi <- summary(fit) ),"try-error")) {
    sufic = fit$coefficients
    coefnames = names(sufic)
    (out = as.numeric(t(sufic)))
    if (coefnames[1] == "(Intercept)")
      (coefnames = c("intercept", paste0("slope",1:(length(out)-1))) )

    names(out) = coefnames

  } else {

    # sufi = summary(fit)
    sufic = sufi$coefficients
    out = as.numeric(t(sufic))

    intercept = attr(sufi$terms,"intercept")

    ## we may want keep forcing through origin, hence no intercept @version 2018-11-29 uncommented 
   # if (nrow(sufic)<2) {
   #  dt = data.table(variable = NA_character_, value = NA_real_)
   #  } else 
    
    ## quantreg can only give a Named num for coefficients @version 2018-11-29
    if (! is.null(ncol(sufic) ) ) { 
    
      if (ncol(sufic) == 4 & intercept == 1){
        # now checking if indeed the t-test is reported
      intnames = c("intercept", "intercept_sd" ,"intercept_ttest", "intercept_pvalue")
      nvars = 1 :(nrow(sufic)-1)
      # nvars = 1:2
      slopenames = c()
        for (nvar in nvars) {
        slopenames = c(slopenames,  paste0("slope",nvar),
          paste0("slope",nvar,"_sd"),
          paste0("slope",nvar,"_ttest"),
          paste0("slope",nvar,"_pvalue"))
          }
        } else if (ncol(sufic) == 4 & intercept == 0){
            intnames = c()
            nvars = 1 :(nrow(sufic))
            # nvars = 1:2
            slopenames = c()
                for (nvar in nvars) {
                slopenames = c(slopenames,  paste0("slope",nvar),
                paste0("slope",nvar,"_sd"),
                paste0("slope",nvar,"_ttest"),
                paste0("slope",nvar,"_pvalue"))
                }
        } else {
        # now for the case of quantreg rq which reports upper and lower bounds
        cna = colnames(sufic)
        cna = sub(pattern = " ", replacement = "", cna)
        cnaend = cna[2:ncol(sufic)]
        (intnames = c("intercept", paste0("intercept_", cnaend)))
        nvars = 1 :(nrow(sufic)-1)
        slopenames = c()
          for (nvar in nvars) {
           slopenames = c(slopenames,  paste0("slope",nvar),
            sapply(cnaend, function(cn) paste0("slope",nvar,"_",cn)) )
          }
          # slopenames
        }
    } else { # capturing the case of coeficients as a named numeric @version 2018-11-29
      intnames = c()
      slopenames =  paste0("slope",(1:length(sufic)) )
    }
      # slopenames
    ## @20160713 error in regression
      #  names(out) = c(intnames, slopenames)
    # Error in names(out) = c(intnames, slopenames) :
    #   'names' attribute [12] must be the same length as the vector [4]

    if (inherits(try(names(out) <- c(intnames, slopenames)),"try-error")) {
        print(paste("ERROR: names of summary object do not match",out, intnames, slopenames))
        print(sufic)
    #     next
      }
    out = c(out,R2adj = sufi$adj.r.squared, n = length(residuals(fit)))
  }

  # out = c(out,R2 = 1 - var(fit$residuals,na.rm = TRUE) / var(fit$y, na.rm = TRUE) )  ### fit$y is only returned by rq() not by lm()
  out = c(out,R2 = 1 - var(fit$residuals,na.rm = TRUE) / var(fit$model[,1], na.rm = TRUE) )

  dt = data.table(statistic = names(out), value = as.numeric(out))
 return(dt)
  }
}


mlm.output.statlong.call = function(mula, data, ...) {
#' Call a linear regression model and reorganize its output 
#' 
#' sucessfull error handling for errors on function calls when groups are empty
#' 
#' @param mula a formula for a lm() regression provided as character using the column names of the data
#' @param data the data table with the column names
#' @param ... additional arguments to call lm() such as weights and other lm arguments 
#' @author Maik Renner, mrenner [at] bgc-jena.mpg.de
#' @return a data.table with columns statistic and value
#' @examples
#' DT = as.data.table(mtcars)
#' DT[, mlm.output.statlong.call(mula = "mpg ~ drat + gear", data = .SD) , by = list(am) ]
#' # impose an error with a nonexisting column name
#' DT[, mlm.output.statlong.call(mula = "mpg ~ drat + gear + bogus", data = .SD) , by = list(am) ]
#' # use ... to weight the regression
#' DT[, mlm.output.statlong.call(mula = "mpg ~ drat + gear", data = .SD, weight = cyl) , by = list(am) ]

#' @version 0.10 2016-05-09  missing value handling in case everything is NA
#' @version 0.11 2017-09-11 check that data is sufficiently long
#' @version 0.12 2017-09-11 an error occurend when all data was the same, then lm worked but yielded an NA for the slopes , thats why another else if statement was used
#' @version 0.13 2018-02-02 change variable to statistic to handle melt operations
#' @version 0.14 2018-03-05 set try to silent to avoid printout on bash
  if (is.null(data) |  nrow(data) < 3 |  inherits(try(ans<-lm(formula = mula, data, ...),silent = TRUE),"try-error")) {
#     print(paste(mula, "combi no data", collapse = " "))
### some non null output for groups is still required !
# otherwise j doesn't evaluate to the same number of columns for each group
    data.table(statistic = NA_character_, value = NA_real_)

  } else if (any(is.na(coef(ans)))) {
    data.table(statistic = NA_character_, value = NA_real_)
  } else {
    mlm.output.statlong(ans)
    }
}


mlm.output.statlong.call.rq = function(mula, data, ...) {
  #' Call a linear QUANTILE regression model and reorganize its output 
  #' 
  #' sucessfull error handling for errors on function calls when groups are empty
  #' 
  #' @param mula a formula for a lm() regression provided as character using the column names of the data
  #' @param data the data table with the column names
  #' @param ... additional arguments to call rq() in particular the quantile must be set e.g. tau = .90 
  #' @author Maik Renner, mrenner [at] bgc-jena.mpg.de
  #' @return a data.table with columns statistic and value
  #' @examples
  #' DT = as.data.table(mtcars)
  #' DT[, mlm.output.statlong.call(mula = "mpg ~ drat + gear", data = .SD) , by = list(am) ]
  #' # impose an error with a nonexisting column name
  #' DT[, mlm.output.statlong.call(mula = "mpg ~ drat + gear + bogus", data = .SD) , by = list(am) ]
  #' # use ... to weight the regression
  #' DT[, mlm.output.statlong.call(mula = "mpg ~ drat + gear", data = .SD, weight = cyl) , by = list(am) ]
  
  #' @version 0.14 2018-03-05 set try to silent to avoid printout on bash
  #' @version 0.15 2018-11-20 based on mlm.output.statlong.call with little adaptions to run rq()
  if (is.null(data) |  nrow(data) < 3 |  inherits(try(ans<-rq(formula = mula, data, ...),silent = TRUE),"try-error")) {
    #     print(paste(mula, "combi no data", collapse = " "))
    ### some non null output for groups is still required !
    # otherwise j doesn't evaluate to the same number of columns for each group
    data.table(statistic = NA_character_, value = NA_real_)
    
  } else if (any(is.na(coef(ans)))) {
    data.table(statistic = NA_character_, value = NA_real_)
  } else {
    mrq.output.statlong(ans)
  }
}



# DT = as.data.table(mtcars)
# DT[, mlm.output.statlong.call(mula = "mpg ~ drat + gear", data = .SD) , by = list(am) ]
# # impose an error with a nonexisting column name
# DT[, mlm.output.statlong.call(mula = "mpg ~ drat + gear + bogus", data = .SD) , by = list(am) ]
# # use ... to weight the regression
# DT[, mlm.output.statlong.call(mula = "mpg ~ drat + gear", data = .SD, weight = cyl) , by = list(am) ]


mrq.output.statlong = function(fit) {
  #' Reorganize the output of a quantile regression fit into a statistic | value format 
  #' 
  #' simple function which takes a lm model output to create a named vector intended for use with data.table and group by
  #' @param fit object returned from a linear model quantreg::rq() are supported
  #' @return a data.table with columns statistic and value
  #' @author Maik Renner, mrenner [at] bgc-jena.mpg.de
  #' @examples
  #' # output of regression model is melted to two columns
  #' library(quantreg)
  #' data(engel)
  #' fit = rq(foodexp~income,tau=0.9,data=engel)
  #' mrq.output.statlong(fit)

  #' @version 0.01 2018-12-04 based on mlm.output.statlong()  
  #'                removing the R2 and adding the pseudo R1 
  #' 
  #' 
  
  if (! is.null(fit) ) {
    if (inherits(try(sufi <- summary(fit) ),"try-error")) {
      sufic = fit$coefficients
      coefnames = names(sufic)
      (out = as.numeric(t(sufic)))
      if (coefnames[1] == "(Intercept)")
        (coefnames = c("intercept", paste0("slope",1:(length(out)-1))) )
      
      names(out) = coefnames
      
    } else {
      
      # sufi = summary(fit)
      sufic = sufi$coefficients
      out = as.numeric(t(sufic))
      
      intercept = attr(sufi$terms,"intercept")
      
      ## we may want keep forcing through origin, hence no intercept @version 2018-11-29 uncommented 
      # if (nrow(sufic)<2) {
      #  dt = data.table(variable = NA_character_, value = NA_real_)
      #  } else 
      
      ## quantreg can only give a Named num for coefficients @version 2018-11-29
      if (! is.null(ncol(sufic) ) ) { 
        
        if (ncol(sufic) == 4 & intercept == 1){
          # now checking if indeed the t-test is reported
          intnames = c("intercept", "intercept_sd" ,"intercept_ttest", "intercept_pvalue")
          nvars = 1 :(nrow(sufic)-1)
          # nvars = 1:2
          slopenames = c()
          for (nvar in nvars) {
            slopenames = c(slopenames,  paste0("slope",nvar),
                           paste0("slope",nvar,"_sd"),
                           paste0("slope",nvar,"_ttest"),
                           paste0("slope",nvar,"_pvalue"))
          }
        } else if (ncol(sufic) == 4 & intercept == 0){
          intnames = c()
          nvars = 1 :(nrow(sufic))
          # nvars = 1:2
          slopenames = c()
          for (nvar in nvars) {
            slopenames = c(slopenames,  paste0("slope",nvar),
                           paste0("slope",nvar,"_sd"),
                           paste0("slope",nvar,"_ttest"),
                           paste0("slope",nvar,"_pvalue"))
          }
        } else {
          # now for the case of quantreg rq which reports upper and lower bounds
          cna = colnames(sufic)
          cna = sub(pattern = " ", replacement = "", cna)
          cnaend = cna[2:ncol(sufic)]
          (intnames = c("intercept", paste0("intercept_", cnaend)))
          nvars = 1 :(nrow(sufic)-1)
          slopenames = c()
          for (nvar in nvars) {
            slopenames = c(slopenames,  paste0("slope",nvar),
                           sapply(cnaend, function(cn) paste0("slope",nvar,"_",cn)) )
          }
          # slopenames
        }
      } else { # capturing the case of coeficients as a named numeric @version 2018-11-29
        intnames = c()
        slopenames =  paste0("slope",(1:length(sufic)) )
      }
      # slopenames
      ## @20160713 error in regression
      #  names(out) = c(intnames, slopenames)
      # Error in names(out) = c(intnames, slopenames) :
      #   'names' attribute [12] must be the same length as the vector [4]
      
      if (inherits(try(names(out) <- c(intnames, slopenames)),"try-error")) {
        print(paste("ERROR: names of summary object do not match",out, intnames, slopenames))
        print(sufic)
        #     next
      }
      out = c(out, n = length(residuals(fit)))
    }
    #The pseudo-R^2 measure suggested by Koenker and Machado's 1999 JASA paper
    # measures goodness of fit by comparing the sum of weighted deviations
    # for the model of interest with the same sum from a model in which only the intercept appears.
    
    R1 = 1 - fit$rho / rq(fit$y ~ 1)$rho
    out = c(out,R1 = R1 )
    
    dt = data.table(statistic = names(out), value = as.numeric(out))
    return(dt)
  }
}
