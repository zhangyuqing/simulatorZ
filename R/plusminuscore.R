### filename: plusminuscore.r
### Title: Backend for plusMinusSurv.r
###
### Author: Levi Waldron
### email: <lwaldron.research@gmail.com>
### date of creation: Mar 21, 2011
###
##### Brief description:
###
###  modeltype="plusminus": returns coefficients of a simple
###  regression model with coefficients of +/- 1, optionally scaled by
###  the standard deviation of the variables.
###
###  modeltype="compoundcovariate": also supports the Compound Covariate method,
###  where coefficients
###  are equal to the univariate Cox regression coefficient.
###
###  modeltype="voting": Each feature receives one vote for positive or negative 
###	risk.
###  
###  modeltype="positiveriskvoting": Each feature can only vote for positive 
###	risk (poor prognosis) or abstain.
###  
###  modeltype="negativeriskvoting": Each feature can only vote for negative 
###	risk (good prognosis) or abstain.
###  
###  **************************************************************************
###  

plusMinus <- structure(function
### function for plusMinus algorithm
( X, 
### gene expression matrix
y, 
### response variables
lambda=NULL, tuningpar="nfeatures",
                       standardize=FALSE, directionality="posneg",
                       ties.method="average",
                       votingthresholdquantile=0.5, modeltype="plusminus")
{
  supported.modeltypes <- c("plusminus", "compoundcovariate", "tscore", 
                            "voting", "positiveriskvoting", "negativeriskvoting")
  if(!modeltype %in% supported.modeltypes)
    stop(paste("provided modeltype:", modeltype, "is not one of the 
               currently supported modeltypes:", 
               supported.modeltypes))
  unitest.output <- rowCoxTests( X=t(X), y=y )
  if(identical(modeltype, "compoundcovariate")){
    cc <- unitest.output$coef
    standardize <- FALSE
  }else{
    cc <- sign(unitest.output$coef)
  }
  if(identical(modeltype, "tscore") & directionality != "posneg")
    stop("modeltype=tscore is compatible only with directionality=posneg.")
  names(cc) <- colnames(X)
  if(identical(directionality,"pos")){
    unitest.output$p.value[cc < 0] <- NA
    cc[cc < 0] <- 0
  }else if(identical(directionality,"neg")){
    unitest.output$p.value[cc > 0] <- NA
    cc[cc > 0] <- 0
  }else if(!identical(directionality,"posneg")){
    stop("directionality must be pos, neg, or posneg")
  }
  if(identical(tuningpar,"nfeatures")){
    unitest.order <- rank(unitest.output$p.value,ties.method=ties.method)
    keep.features <- unitest.order <= lambda
  }else if(identical(tuningpar,"pval")){
    keep.features <- unitest.output$p.value <= lambda
  }else{
    stop("tuningpar must be nfeatures or pval")
  }
  cc[!keep.features] <- 0
  if(standardize){
    cc <- cc / apply(X,2,sd)
  }
  ##set thresholds for voting schemes - each feature will cast its
  ##vote according to whether it is above or below this threshold.
  if (modeltype %in% c("voting", "positiveriskvoting", "negativeriskvoting")){
    ##Note that these three voting schemes are differentiated only in their
    ###predict methods.
    if (standardize)
      warning("Standardize would normally be FALSE in a voting scheme; 
                            setting standardize=TRUE will weight the votes 
                            by the Cox  coefficient of each feature.")
    votingthresholds <- apply(X, 2, quantile, probs=votingthresholdquantile)
    ##If cc is zero, set threshold to zero:
    votingthresholds[abs(cc) < .Machine$double.eps ^ 0.5] <- 0
    ret.obj <- new("ModelLinear", coefficients = cc, modeltype = modeltype, 
                   votingthresholds=votingthresholds)
  }else{
    ##scale coefficients for numerical stability simply by
    ##dividing by the number of non-zero coefficients.
    cc <- cc / sum(abs(cc) > 0)
    #ret.obj <- new("ModelLinear", coefficients = cc, modeltype = modeltype)
    ret.obj <- cc
  }
  return(ret.obj)
  ### returns regression coefficients 
},ex=function(){
  set.seed(8)
  library(curatedOvarianData)
  data( E.MTAB.386_eset )
  eset <- E.MTAB.386_eset[1:100, 1:30]
  
  X <- t(exprs(eset))
  
  time <- eset$days_to_death
  cens <- sample(0:1, 30, replace=TRUE)
  y <- Surv(time, cens)
  
  beta <- plusMinus(X, y)
  beta
})