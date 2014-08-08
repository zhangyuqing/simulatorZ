simData <- function
### Input is a list of ExpressionSets of the original sets, 
### output is another list of ExpressionSets of the simulated sets
### Contains the two steps of non-parametric bootstrap 
(esets,
 ### original data sets
 balance.variables,
 ### Later, balance.variables will be a vector of variable names that should be balanced in the simulation.
 n.samples,
 ### how many samples should be generated
 type
 ### original or simulated
){  
  #### Step.1 Drawing sets
  if(type != "original" && type != "simulated")
    stop("Wrong type.")
  else if(type == "original")
    setsID <- 1:length(esets)
  else{
    prob.set <- rep((1/length(esets)), times=length(esets))
    setsID <- sample(1:length(esets), prob=prob.set, replace=TRUE)  # labels of data sets
  }  
  print(setsID) 
  
  #### Step.2 Calculate the joint probability distribution  
  ## rbinding all covariates data  
  covariates.list <- lapply(esets, function(eset){
    if(length(balance.variables) == 1){
      return(as.character(pData(eset)[, balance.variables]))
    }
    else{
      return(as.character(do.call(paste, pData(eset)[, balance.variables])))
    }
  })
  covariate_all <- do.call(c, covariates.list)
  prob_desired <- table(covariate_all) / sum(table(covariate_all))  
  probs.list <- lapply(covariates.list, function(covariates){
    prob_real <- table(covariates) / sum(table(covariates))
    prob_desired_matched <- prob_desired[match(names(prob_real), names(prob_desired))]
    if(identical(names(prob_desired_matched), names(prob_real)))
      prob <- prob_desired_matched / prob_real
    sample_probs <- prob[match(covariates, names(prob))]
    sample_probs <- sample_probs / sum(sample_probs)
    return(sample_probs)
  })  
  #### Step. 3 Drawing patients
  samplesets <- sampleind <- list()
  if(!is.null(balance.variables)){      
    print(paste("covariate: ", balance.variables, sep=""))
    for(i in 1:length(esets)){
      sampleind[[i]] <- sample(1:length(sampleNames(esets[[setsID[i]]])), n.samples, replace=TRUE, prob=probs.list[[setsID[i]]])
      samplesets[[i]] <- esets[[setsID[i]]][, sampleind[[i]]]
    }        
  }    
  else {
    print("covariate: NULL")
    for(i in 1:length(esets)){
      sampleind[[i]] <- sample(1:length(sampleNames(esets[[setsID[i]]])), n.samples, replace=TRUE)
      samplesets[[i]] <- esets[[setsID[i]]][, sampleind[[i]]]
    }
  }
  
  names(samplesets) <- setsID
  ## We still need to remember which studies these samples are from, because we will need to 
  ## generate survival time later, and beta from CoxBoost should be relevant to original data sets 
  ## End Step.2
  res <- list(esets=samplesets, indices=sampleind, setsID=setsID)
  return(res)
  ### returns a list of simulated ExpressionSets, with names indicating its origin.
}