simTime <- structure(function
### simTime is a function to perform the parametric-bootstrap step, where we use the true coefficients
### and cumulative hazard to simulate survival and censoring.
(simmodels,
 ### a list in the form of the return value of simData()
 ### which consists of three lists:
 ### esets: a list of ExpressionSets
 ### setsID: a list of set labels indicating which original set the simulated one is from
 ### indices: a list of patient labels to tell which patient in the original set is drawn
 result
 ### a list in the form of return of getTrueModel()
 ### which consists of five lists: 
 ### Beta: a list of coefficients obtained by 
 ### grid: timeline grid corresponding to hazard estimations censH and survH
 ### survH: cumulative hazard for survival times distribution
 ### censH: cumulative hazard for censoring times distribution
 ### lp: true linear predictors
){  
  y.vars <- list()
  for(i in 1:length(simmodels$esets)){
    print(i)
    setID <- simmodels$setsID[i]
    time <- simmodels$y.vars[[setID]][, 1]
    time <- as.numeric(as.character(time))
    status <-  simmodels$y.vars[[setID]][, 2] 
    status <- as.numeric(as.character(status))
    ith_beta <- result$beta[[setID]]
    grid <- result$grid[[setID]]
    survH <- result$survH[[setID]]
    censH <- result$censH[[setID]]
    lp <- result$lp[[setID]][simmodels$indices[[i]]]
    ## get the inverse of the Nelsin-Aalen estimator function
    ## the ith resampled set is from names(esets)[i] originally(refer to simData())
    TIME <- CENS <- c()
    newTIME <- newCENS <- c()
    for(j in 1:length(sampleNames(simmodels$esets[[i]]))) {
      ## generate the survival time 
      u1 <- runif(1, min = 0, max = 1)
      z1 <- (-log(u1, base=exp(1))) * (exp(-lp[j])) 
      TIME[j] <- grid[which.min(abs(survH - z1))]
      if (TIME[j] >= max(time[which(status == 1)], na.rm=TRUE))  TIME[j] <- 1e+08      
      ## generate the censoring time
      u2 <- runif(1, min = 0, max = 1)
      z2 <- -log(u2, base=exp(1))
      CENS[j] <- min(grid[which.min(abs(censH - z2))], max(time[which(status == 0)]))
      if (CENS[j] >= max(time[which(status == 0)], na.rm=TRUE))  CENS[j] <- 1e+08            
      newTIME[j] <- min(TIME[j], CENS[j])
      newCENS[j] <- as.numeric(TIME[j] < CENS[j])
    }    
    y.vars[[i]] <- Surv(newTIME, newCENS)
  }
  simmodels$y.vars <- y.vars
  return(simmodels)
  ### survival time is saved in phenodata, here the function still returns the ExpressionSets
},ex=function(){
  library(curatedOvarianData)
  data( E.MTAB.386_eset )
  eset1 <- E.MTAB.386_eset[1:10, 1:5]
  eset2 <- E.MTAB.386_eset[1:10, 6:10]
  eset3 <- E.MTAB.386_eset[1:10, 11:15]  
  esets <- list(eset1, eset2, eset3) 
  
  time1 <- eset1$days_to_death
  cens1 <- c(0, 0, 0, 1, 1)
  y1 <- Surv(time1, cens1)
  time2 <- eset2$days_to_death
  cens2 <- c(1, 1, 0, 0, 0)
  y2 <- Surv(time2, cens2)
  time3 <- eset3$days_to_death
  cens3 <- c(1, 0, 0, 0, 1)
  y3 <- Surv(time3, cens3)
  y.vars <- list(y1, y2, y3)
  
  # To perform both parametric and non-parametric bootstrap, you can call simBootstrap()
  # or, you can divide the steps into:
  res <- getTrueModel(esets, y.vars, 100)
  simmodels <- simData(esets=esets, y.vars=y.vars, n.samples=10)
  
  # Then, use this function
  simmodels <- simTime(simmodels=simmodels, result=res) 
  
  # it also supports performing only the parametrc bootstrap step on a list of expressionsets
  # but you need to construct the parameter by scratch
  res <- getTrueModel(esets, y.vars, 100)
  setsID <- 1:length(esets)
  indices <- list()
  for(i in setsID){
    indices[[i]] <- 1:length(sampleNames(esets[[i]])) 
  }
  simmodels <- list(esets=esets, y.vars=y.vars, indices=indices, setsID=setsID)
  
  new.simmodels <- simTime(simmodels=simmodels, result=res)  
})