simData <- function
### Input is a list of ExpressionSets of the original sets, 
### output is another list of ExpressionSets of the simulated sets
### Contains the two steps of non-parametric bootstrap 
(esets,
 ### original data sets
 balance.variables
 ### Later, balance.variables will be a vector of variable names that should be balanced in the simulation.
){  
  #### Step.1 Drawing sets
  prob.set <- rep((1/length(esets)), times=length(esets))
  setsID <- sample(1:length(esets), prob=prob.set, replace=TRUE)  # labels of data sets
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
  samplesets <- list()
  if(!is.null(balance.variables)){      
    print(paste("covariate: ", balance.variables, sep=""))
    for(i in 1:length(esets)){
      sampleind <- sample(1:length(sampleNames(esets[[setsID[i]]])), 150, replace=TRUE, prob=probs.list[[setsID[i]]])
      samplesets[[i]] <- esets[[setsID[i]]][, sampleind]
    }        
  }    
  else {
    print("covariate: NULL")
    for(i in 1:length(esets)){
      sampleind <- sample(1:length(sampleNames(esets[[setsID[i]]])), 150, replace=TRUE)
      samplesets[[i]] <- esets[[setsID[i]]][, sampleind]
    }
  }
  
  names(samplesets) <- setsID
  ## We still need to remember which studies these samples are from, because we will need to 
  ## generate survival time later, and beta from CoxBoost should be relevant to original data sets 
  ## End Step.2
  return(samplesets)
  ### returns a list of simulated ExpressionSets, with names indicating its origin.
}