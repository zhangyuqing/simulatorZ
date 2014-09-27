getTrueModel <- structure(function
### The parametric bootstrap simulation depends on the true model of original sets.
### This function is to generate useful values from the true models for further analysis.  
### We fit CoxBoost to the original sets and use the coefficients to simulate
### the survival and censoring time. grid, survH, censH, which are useful for this purpose.
### grid=grid corresponding to hazard estimations censH and survH
### survH=cumulative hazard for survival times distribution
### censH=cumulative hazard for censoring times distribution
(obj,
 ### a list of ExpressionSets, matrix or SummarizedExperiment
 y.vars,
 ### a list of response variables, Surv, matrix or data.frame object
 parstep 
 ### number of steps in CoxBoost
){
  Beta <- grid <- survH <- censH <- result <- lp <- list()
  for(i in 1:length(obj)){
    print(i)
    ### PART 1: get beta
    time <- y.vars[[i]][, 1]
    time <- as.numeric(as.character(time))
    status <- y.vars[[i]][, 2]
    status <- as.numeric(as.character(status))
    
    if(class(obj[[1]])=="ExpressionSet"){
      X <- t(exprs(obj[[i]]))
    }      
    else if(class(obj[[1]])=="matrix"){
      X <- t(obj[[i]])
    }  
    else if(class(obj[[1]])=="SummarizedExperiment"){
      X <- t(assay(obj[[i]]))
    }
      
      
    cbfit <- CoxBoost(time=time, status=status, x=X, 
                      stepno=parstep, standardize=FALSE)
    
    Beta[[i]] <- coef(cbfit)    
    ### PART 2: get grid, suvH, censH
    ### Survival time ==> dmfs.cens=1, Censoring time ==> dmfs.cens=0
    ### get survH, censH, grid    
    
    if(class(obj[[1]])=="ExpressionSet")
      X <- exprs(obj[[i]])
    else if(class(obj[[1]])=="matrix")
      X <- obj[[i]]
    else if(class(obj[[1]])=="SummarizedExperiment")
      X <- assay(obj[[i]])
    
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
  data(GSE17260_eset)
  data(E.MTAB.386_eset)
  data(GSE14764_eset)
  esets <- list(GSE17260=GSE17260_eset, E.MTAB.386=E.MTAB.386_eset, GSE14764=GSE14764_eset)
  esets.list <- lapply(esets, function(eset){
    return(eset[1:500, 1:20])
  })
  
  ## simulate on multiple ExpressionSets
  set.seed(8) 
  
  y.list <- lapply(esets.list, function(eset){
    time <- eset$days_to_death
    cens.chr <- eset$vital_status
    cens <- c()
    for(i in 1:length(cens.chr)){
      if(cens.chr[i] == "living") cens[i] <- 1
      else cens[i] <- 0
    }
    y <- Surv(time, cens)
    return(y)
  })
     
  res1 <- getTrueModel(esets.list, y.list, 100)
  ## Get true model from one set
  res2 <- getTrueModel(esets.list[1], y.list[1], 100)
  names(res2)
  res2$lp
  ## note that y.list[1] cannot be replaced by y.list[[1]]
   
})