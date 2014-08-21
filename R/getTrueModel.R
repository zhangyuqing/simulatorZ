getTrueModel <- structure(function
### The parametric bootstrap simulation depends on the true model of original sets.
### This function is to generate useful values from the true models for further analysis.  
### We fit CoxBoost to the original sets and use the coefficients to simulate
### the survival and censoring time. grid, survH, censH, which are useful for this purpose.
### grid=grid corresponding to hazard estimations censH and survH
### survH=cumulative hazard for survival times distribution
### censH=cumulative hazard for censoring times distribution
(esets,
 ### a list of ExpressionSets
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
    cbfit <- CoxBoost(time=time, status=status, x=t(exprs(esets[[i]])), 
                      stepno=parstep, standardize=FALSE)
    
    Beta[[i]] <- coef(cbfit)    
    ### PART 2: get grid, suvH, censH
    ### Survival time ==> dmfs.cens=1, Censoring time ==> dmfs.cens=0
    ### get survH, censH, grid    
    X <- exprs(esets[[i]])
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
  data( E.MTAB.386_eset )
  eset1 <- E.MTAB.386_eset[, 1:5]
  eset2 <- E.MTAB.386_eset[, 6:10]
  eset3 <- E.MTAB.386_eset[, 11:15]
  
  ## simulate on multiple ExpressionSets
  esets.list <- list(eset1, eset2, eset3) 
  
  y.list <- list()
  for(i in 1:length(esets.list)){
    time <- esets.list[[i]]$days_to_death
    cens <- sample(0:1, 5, replace=TRUE)
    y.list[[i]] <- Surv(time, cens)
  }  
  
  res1 <- getTrueModel(esets.list, y.list, 100)
  res2 <- getTrueModel(list(eset1), y.list[1], 100)
  ## note that y.list[1] cannot be replaced by y.list[[1]]
  
})