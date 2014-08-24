getTrueModel <- structure(function
### The parametric bootstrap simulation depends on the true model of original sets.
### This function is to generate useful values from the true models for further analysis.  
### We fit CoxBoost to the original sets and use the coefficients to simulate
### the survival and censoring time. grid, survH, censH, which are useful for this purpose.
### grid=grid corresponding to hazard estimations censH and survH
### survH=cumulative hazard for survival times distribution
### censH=cumulative hazard for censoring times distribution
(esets,
 ### a list of ExpressionSets, matrix or SummarizedExperiment
 y.vars,
 ### a list of response variables
 parstep 
 ### CoxBoost parameter
 #par_penalty
 ### CoxBoost parameter
){
  Beta <- grid <- survH <- censH <- result <- lp <- list()
  for(i in 1:length(esets)){
    print(i)
    ### PART 1: get beta
    time <- y.vars[[i]][, 1]
    time <- as.numeric(as.character(time))
    status <- y.vars[[i]][, 2]
    status <- as.numeric(as.character(status))
    
    if(class(esets[[1]])=="ExpressionSet"){
      X <- t(exprs(esets[[i]]))
    }      
    else if(class(esets[[1]])=="matrix"){
      X <- t(esets[[i]])
    }  
    else if(class(esets[[1]])=="SummarizedExperiment"){
      X <- t(assay(esets[[i]]))
    }
      
      
    cbfit <- CoxBoost(time=time, status=status, x=X, 
                      stepno=parstep, standardize=FALSE)
    
    Beta[[i]] <- coef(cbfit)    
    ### PART 2: get grid, suvH, censH
    ### Survival time ==> dmfs.cens=1, Censoring time ==> dmfs.cens=0
    ### get survH, censH, grid    
    
    if(class(esets[[1]])=="ExpressionSet")
      X <- exprs(esets[[i]])
    else if(class(esets[[1]])=="matrix")
      X <- esets[[i]]
    else if(class(esets[[1]])=="SummarizedExperiment")
      X <- assay(esets[[i]])
    
    lp[[i]] <- as.numeric(Beta[[i]] %*% X)  ## Calculate linear predictor
    grid[[i]] <- seq(0, max(time), by = 1)
    
    survH[[i]] <- basehaz.gbm(t=time, delta=status, f.x=lp[[i]], 
                              t.eval=grid[[i]], smooth=TRUE, cumulative=TRUE) 
    
    inverse_status <- (-1) * status + 1
    censH[[i]] <- basehaz.gbm(t=time, delta=inverse_status, f.x=rep(0, length(time)), 
                              t.eval=grid[[i]], smooth=TRUE, cumulative=TRUE)
    
  }
  result <- list(beta=Beta, grid=grid, survH=survH, censH=censH, lp=lp)  
  return(result)
  ### returns a list of values:
  ### beta: True coefficients obtained by fitting CoxBoost to the original ExpressionSets
  ### grid: timeline grid corresponding to hazard estimations censH and survH
  ### survH: cumulative hazard for survival times distribution
  ### censH: cumulative hazard for censoring times distribution
  ### lp: true linear predictors 
},ex=function(){
  library(curatedOvarianData)
  library(GenomicRanges)
  data( E.MTAB.386_eset )
  eset1 <- E.MTAB.386_eset[, 1:5]
  eset2 <- E.MTAB.386_eset[, 6:10]
  eset3 <- E.MTAB.386_eset[, 11:15]
  
  ## simulate on multiple ExpressionSets
  esets.list <- list(eset1, eset2, eset3) 
  
  time1 <- eset1$days_to_death
  cens1 <- c(0, 0, 0, 1, 1)
  y1 <- Surv(time1, cens1)
  time2 <- eset2$days_to_death
  cens2 <- c(1, 1, 0, 0, 0)
  y2 <- Surv(time2, cens2)
  time3 <- eset3$days_to_death
  cens3 <- c(1, 0, 0, 0, 1)
  y3 <- Surv(time3, cens3)
  y.list<- list(y1, y2, y3) 
     
  res1 <- getTrueModel(esets.list, y.list, 100)
  ## Get true model from one set
  res2 <- getTrueModel(list(eset1), y.list[1], 100)
  names(res2)
  res2$lp
  ## note that y.list[1] cannot be replaced by y.list[[1]]
  
  ## Support matrices
  X.list <- lapply(esets.list, function(eset){
    return(exprs(eset))
  })
  res3 <- getTrueModel(X.list, y.list, 100)
  
})