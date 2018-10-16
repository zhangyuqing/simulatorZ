simTime <- structure(function
### simTime is a function to perform the parametric-bootstrap step, where we use the true coefficients
### and cumulative hazard to simulate survival and censoring.
(simmodels,
 ### a list in the form of the return value of simData()
 ### which consists of three lists:
 ### obj: a list of ExpressionSets, matrices or RangedSummarizedExperiments
 ### setsID: a list of set labels indicating which original set the simulated one is from
 ### indices: a list of patient labels to tell which patient in the original set is drawn
 original.yvars,
 ### response variable in the order of original sets(without sampling)
 result
 ### a list in the form of return of getTrueModel()
 ### which consists of five lists: 
 ### Beta: a list of coefficients obtained by 
 ### grid: timeline grid corresponding to hazard estimations censH and survH
 ### survH: cumulative hazard for survival times distribution
 ### censH: cumulative hazard for censoring times distribution
 ### lp: true linear predictors
){  
  new.y.vars <- list()
  for(i in seq_along(simmodels$obj)){
    print(i)
    setID <- simmodels$setsID[i]
    time <- original.yvars[[setID]][, 1]
    time <- as.numeric(as.character(time))
    status <-  original.yvars[[setID]][, 2] 
    status <- as.numeric(as.character(status))
    ith_beta <- result$beta[[setID]]
    grid <- result$grid[[setID]]
    survH <- result$survH[[setID]]
    censH <- result$censH[[setID]]
    lp <- result$lp[[setID]][simmodels$indices[[i]]] ## bootstrap of true linear scores
    ## get the inverse of the Nelsin-Aalen estimator function
    ## the ith resampled set is from names(esets)[i] originally(refer to simData())
    TIME <- CENS <- c()
    newTIME <- newCENS <- c()
    
    num.sam <- ncol(getMatrix(simmodels$obj[[i]])) 
    for(j in seq_len(num.sam)) {
      ## generate the survival time 
      u1 <- runif(1, min = 0, max = 1)
      z1 <- -log(u1, base=exp(1)) * exp(-lp[j])
      TIME[j] <- grid[which.min(abs(survH - z1))]
      if (TIME[j] >= max(time[which(status == 1)], na.rm=TRUE))  TIME[j] <- 1e+08      
      ## generate the censoring time
      u2 <- runif(1, min = 0, max = 1)
      z2 <- -log(u2, base=exp(1)) ## * exp(-lp) and lp=0
      CENS[j] <- min(grid[which.min(abs(censH - z2))], max(time[which(status == 0)]))
      #if (CENS[j] >= max(time[which(status == 0)], na.rm=TRUE))  CENS[j] <- 1e+08            
      newTIME[j] <- min(TIME[j], CENS[j])
      newCENS[j] <- as.numeric(TIME[j] < CENS[j])
    }    
    new.y.vars[[i]] <- Surv(newTIME, newCENS)
  }
  simmodels$y.vars <- new.y.vars
  return(simmodels)
  ### survival time is saved in phenodata, here the function still returns the ExpressionSets
},ex=function(){
  library(curatedOvarianData)
  data(E.MTAB.386_eset)
  data(GSE14764_eset)
  esets.list <- list(E.MTAB.386=E.MTAB.386_eset[1:100, 1:20], GSE14764=GSE14764_eset[1:100, 1:20])
  rm(E.MTAB.386_eset, GSE14764_eset)
  
  ## simulate on multiple ExpressionSets
  set.seed(8) 
  
  y.list <- lapply(esets.list, function(eset){
    time <- eset$days_to_death
    cens.chr <- eset$vital_status
    cens <- rep(0, length(cens.chr))
    cens[cens.chr=="living"] <- 1
    return(Surv(time, cens))
  })
  
  # To perform both parametric and non-parametric bootstrap, you can call simBootstrap()
  # or, you can divide the steps into:
  res <- getTrueModel(esets.list, y.list, 100)
  simmodels <- simData(obj=esets.list, y.vars=y.list, n.samples=10)
  
  # Then, use this function
  simmodels <- simTime(simmodels=simmodels, original.yvars=y.list, result=res) 
  
  # it also supports performing only the parametrc bootstrap step on a list of expressionsets
  # but you need to construct the parameter by scratch
  res <- getTrueModel(esets.list, y.list, 100)
  setsID <- seq_along(esets.list)
  indices <- list()
  for(i in setsID){
    indices[[i]] <- seq_along(sampleNames(esets.list[[i]])) 
  }
  simmodels <- list(obj=esets.list, y.vars=y.list, indices=indices, setsID=setsID)
  
  new.simmodels <- simTime(simmodels=simmodels, original.yvars=y.list, result=res)  
})
