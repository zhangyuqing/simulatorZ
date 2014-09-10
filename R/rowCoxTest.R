# Filename: filter.r
# Title: Gene selection (filter) methods.
#
# Author: Christoph Bernau
# Email: <bernau@ibe.med.uni-muenchen.de>
# Date of creation: 28.8.2012
#
# Brief description:
#   Returns an object of class 'GeneSel'.
#
# Further comments and notes:
#   Are usually not called directly by the user, but via
#   'GeneSelection'.
#
###############################################################################

coxmatC<-function
### function to perform fast Cox regression
(X,time,status){
  method <- "efron"
  result<-  .C("coxmat", regmat = as.double(X), ncolmat = 
                 as.integer(ncol(X)), nrowmat = as.integer(nrow(X)), 
               reg = as.double(X[, 1]), zscores = as.double(numeric(ncol(X))),
               coefs = as.double(numeric(ncol(X))), maxiter = as.integer(20), 
               nusedx = as.integer(nrow(X)), nvarx = as.integer(1), 
               time = as.double(time), status = as.integer(status), 
               offset = as.double(numeric(nrow(X))), 
               weights = as.double(numeric(nrow(X)) + 1), 
               strata = as.integer(numeric(nrow(X))), means = double(1), 
               beta = double(1), u = double(1), imat = double(1), 
               loglik = double(2), flag = integer(1), work = 
                 double(2 * nrow(X) + 2 + 3), eps = as.double(1e-09), 
               tol_chol = as.double(.Machine$double.eps^0.75), 
               sctest = as.double(method == "efron"), sctest2 = as.double(1), 
               sctest3 = as.double(1))	
  return(result)	
}


fastCox <- function
### Method for performing fast Cox regression
(X, y, learnind, criterion, ...) {
  
  X <- X[learnind, ]
  time <- y[learnind, 1]
  status <- y[learnind, 2]
  sorted <- order(time)
  time <- time[sorted]
  status <- status[sorted]
  X <- as.matrix(X[sorted, ])
  # compute columnwise coxmodels
  out <- coxmatC(X,time,status)
  # compute p-values
  if (criterion == "pvalue") 
    crit <- (1 - pnorm(abs(out$zscores))) * 2
  if (criterion == "coefficient") 
    crit <- abs(out$coefs)
  
  ### and return a VarSelOut-object
  new("VarSelOut", varsel = crit, criterion = criterion)
}

# equivalent to genefilter::rowttests for the cox model.  This is much faster
# than calling coxph for each row of a ##igh-dimensional matrix.
rowCoxTests <- structure(function
### method for performing Cox regression
(X, 
 ### Gene expression data. The following formats are available:   
 ### matrix Rows correspond to observations, columns to variables.
 ### data.frame Rows correspond to observations, columns to variables.
 ### ExpressionSet rowCoxTests will extract the expressions using exprs().
 y, 
 ### Survival Response, an object of class:
 ### Surv if X is of type data.frame or matrix
 ### character if X is of type ExpressionSet. 
 ### In this case y is the name of the survival 
 ### response in the phenoData of X. If survival 
 ### time and indicator are stored separately 
 ### in the phenoData one can specify a two-element 
 ### character vector the first element representing 
 ### the survival time variable.
 option = c("fast", "slow"), 
 ### "fast" loops over rows in C, "slow" calls coxph 
 ### directly in R. The latter method may be used if 
 ### something goes wrong with the "fast" method.
 ...
 ### currently unused
 ) {
 
  option <- match.arg(option)
  if (identical(option, "fast")) {
    X <- t(as.matrix(X))  #make variables columns
    time <- y[, 1]
    status <- y[, 2]
    sorted <- order(time)
    time <- time[sorted]
    status <- status[sorted]
    X <- X[sorted, ]
    ## method for handling ties (alternative 'breslow')
    method <- "efron"
    ## compute columnwise coxmodels
    out <- coxmatC(X,time,status)
    ## compute p-values and return them
    output <- data.frame(coef = out$coefs, se.coef = out$coefs/out$zscores, 
                         p.value = (1 - 
                                      pnorm(abs(out$zscores))) * 2) 
    rownames(output) <- colnames(X)
  } else if (identical(option, "slow")) {
    output <- t(apply(X, 1, function(xrow) {
      fit <- try(coxph(y ~ xrow))
      if (class(fit) == "try-error") {
        c(NA, NA)
      } else {
        summary(fit)$coefficients[1, c(1, 3, 5)]
      }
    }))
    colnames(output) <- c("coef", "se.coef", "p.value")
    rownames(output) <- rownames(X)
    output <- data.frame(output)
  } else stop("rowCoxTests: option should be fast or slow.")
  return(output)
  ### dataframe with two columns: coef = Cox regression
  ###coefficients, p.value =
  ### Wald Test p-values.  Rows correspond to the rows of X.
},ex=function(){
  #test
  ##regressor-matrix (gene expressions)
  X<-matrix(rnorm(1e6),nrow=10000)
  #seed
  set.seed(123)
  #times
  time<-rnorm(n=ncol(X),mean=100)
  #censoring(1->death)
  status<-rbinom(n=ncol(X),size=1, prob=0.8)
  
  ##survival object
  y<-Surv(time,status)
  
  ## Do 10,000 Cox regressions:
  system.time(output <- rowCoxTests(X=X,y=y, option="fast"))
})