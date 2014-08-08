test_simData(){
  load("ExpressionSets_rm_new.RData")
  origin_esets <- esets
  simmodel <- simData(esets=origin_esets, 
                      balance.variables=c("age","size","grade","node"))  
  new_esets <- simmodel$esets
  covariates.origin <- lapply(origin_esets, function(eset){
    if(length(balance.variables) == 1){
      return(as.character(pData(eset)[, balance.variables]))
    }
    else{
      return(as.character(do.call(paste, pData(eset)[, balance.variables])))
    }
  })
  covariates.originall <- do.call(c, covariates.origin)
  covariates.sim <- lapply(new_esets, function(eset){
    if(length(balance.variables) == 1){
      return(as.character(pData(eset)[, balance.variables]))
    }
    else{
      return(as.character(do.call(paste, pData(eset)[, balance.variables])))
    }
  })
  
  for(i in 1:length(covariates.sim)){
    id <- match(covariates.sim[[i]], covariates.originall)
    freq.sim <- table(covariates.sim[[i]]) / sum(table(covariates.sim[[i]]))
    freq.origin <- table(covariates.originall[id]) / sum(table(covariates.originall[id]))
    checkEqualsNumeric(as.numeric(freq.sim), as.numeric(freq.origin))
  }
}

