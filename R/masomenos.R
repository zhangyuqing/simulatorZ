masomenos <- structure(function
### function for Mas-o-menos algorithm
##references<< Zhao, S., Huttenhower, G. P. C., and Waldron, L. (2013). Mas-o-menos:
## a simple sign averaging method for discrimination in genomic data analysis.
## http://biostats.bepress.com/harvardbiostat/paper158/. Accessed: 2013-10-24.
(X,  
 ### matrix with rows corresponding to subjects and columns to features resp
 y,  
 ### response variable, a data.frame, matrix, or Surv object: c(time, event)
 option="fast", 
 ### whether to use C or R code to fit the marginal Cox models
 ...
 ){

  X <- data.frame(X, row.names=NULL)
  if (!is(y, 'Surv')) y <- Surv(y[, 1], y[, 2])
  alpha <- rowCoxTests(t(X), y, option=option, ...) 
  return(sign(alpha$coef) / nrow(alpha))
  ### return the coefficients
},ex=function(){
  set.seed(8)
  library(curatedOvarianData)
  data( E.MTAB.386_eset )
  eset <- E.MTAB.386_eset[1:100, 1:30]
  
  X <- t(exprs(eset))
  
  time <- eset$days_to_death
  cens <- sample(0:1, 30, replace=TRUE)
  y <- Surv(time, cens)
  
  beta <- masomenos(X, y)
  beta
})
