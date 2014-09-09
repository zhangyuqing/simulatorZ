simTime <- structure(function
### simTime is a function to perform the parametric-bootstrap step, where we use the true coefficients
### and cumulative hazard to simulate survival and censoring.
(simmodels,
 ### a list in the form of the return value of simData()
 ### which consists of three lists:
 ### obj: a list of ExpressionSets, matrics or SummarizedExperiments
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
  for(i in 1:length(simmodels$obj)){
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
    
    if(class(simmodels$obj[[1]])=="ExpressionSet")
      num.sam <- ncol(exprs(simmodels$obj[[i]]))
    else if(class(simmodels$obj[[1]])=="matrix")
      num.sam <- ncol(simmodels$obj[[i]])
    else if(class(simmodels$obj[[1]])=="SummarizedExperiment")
      num.sam <- ncol(assay(simmodels$obj[[i]]))
    
    for(j in 1:num.sam) {
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
  library(GenomicRanges)
  source(system.file("extdata", "patientselection.config",
                     package="curatedOvarianData"))
  source(system.file("extdata", "createEsetList.R", package="curatedOvarianData"))
  esets.list <- lapply(esets, function(eset){
    return(eset[1:500, 1:20])
  })
  esets.list <- esets.list[1:5]
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
  
  # To perform both parametric and non-parametric bootstrap, you can call simBootstrap()
  # or, you can divide the steps into:
  res <- getTrueModel(esets.list, y.list, 100)
  simmodels <- simData(obj=esets.list, y.vars=y.list, n.samples=10)
  
  # Then, use this function
  simmodels <- simTime(simmodels=simmodels, result=res) 
  
  # it also supports performing only the parametrc bootstrap step on a list of expressionsets
  # but you need to construct the parameter by scratch
  res <- getTrueModel(esets.list, y.list, 100)
  setsID <- 1:length(esets.list)
  indices <- list()
  for(i in setsID){
    indices[[i]] <- 1:length(sampleNames(esets.list[[i]])) 
  }
  simmodels <- list(obj=esets.list, y.vars=y.list, indices=indices, setsID=setsID)
  
  new.simmodels <- simTime(simmodels=simmodels, result=res)  
})