generateSurvtime <- function
### Step.3 Generating survival time --- parametric bootstrap
### update the dmfs.time in simulated ExpressionSets
(esets, 
 ### simulated ExpressionSets 
 result
 ### return of getParameter(), a list which consists of Beta, grid, survH and censH
 ### for each original data set
){  
  for(i in 1:length(esets)){
    print(i)
    setID <- as.numeric(names(esets)[i])
    time <- pData(esets[[setID]])[, "dmfs.time"] 
    time <- as.numeric(as.character(time))
    status <- pData(esets[[setID]])[, "dmfs.cens"] 
    status <- as.numeric(as.character(status))
    ith_beta <- as.matrix(result$beta[[setID]])
    grid <- result$grid[[setID]]
    survH <- result$survH[[setID]]
    censH <- result$censH[[setID]]
    lp <- as.numeric(t(exprs(esets[[setID]])) %*% ith_beta)
    ## get the inverse of the Nelsin-Aalen estimator function
    ## the ith resampled set is from names(esets)[i] originally(refer to simData())
    TIME <- CENS <- c()
    newTIME <- newCENS <- c()
    for(j in 1:length(sampleNames(esets[[setID]]))) {
      ## generate the survival time 
      u1 <- runif(1, min = 0, max = 1)
      z1 <- (-log(u1, base=exp(1))) * (exp(-lp[j])) 
      TIME[j] <- grid[which.min(abs(survH - z1))]
      if (TIME[j] >= max(time[which(status == 1)], na.rm=TRUE))  TIME[j] <- 1e+08      
      ## generate the censoring time
      u2 <- runif(1, min = 0, max = 1)
      z2 <- -log(u2, base=exp(1))
      CENS[j] <- grid[which.min(abs(censH - z2))]
      if (CENS[j] >= max(time[which(status == 0)], na.rm=TRUE))  CENS[j] <- 1e+08            
      newTIME[j] <- min(TIME[j], CENS[j])
      newCENS[j] <- as.numeric(TIME[j] < CENS[j])
    }
    pData(esets[[setID]])[, "dmfs.time"] <- as.factor(newTIME)
    #print(as.numeric(as.character(pData(phenoData(esets[[i]]))[, "dmfs.time"])))
    pData(esets[[setID]])[, "dmfs.cens"] <- as.factor(newCENS)
    #print(as.numeric(as.character(pData(phenoData(esets[[i]]))[, "dmfs.cens"])))
  }
  return(esets)
  ### survival time is saved in phenodata, here the function still returns the ExpressionSets
}