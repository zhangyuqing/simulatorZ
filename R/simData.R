simData <- structure(function
### simData() is a function to perform non-parametric bootstrap on a list of ExpressionSets. 
### Input is a list of ExpressionSets(we will call it the original ones),  
### output includes a list of the simulated ExpressionSets, as well as the indices of datasets and patients in the original sets. 
### The function supports: 1.choosing whether to resample from sets by changing the parameter type, and 2. choosing whether to balance over part or all of the 
### covariates by changing parameter balance.variables.
(esets,
 ### a list of ExpressionSets
 n.samples,
 ### an integer indicating how many samples should be resampled from each set
 y.vars=list(),
 ### a list of response variables
 type="two-steps",
 ### "one-step" or "two-steps". If type="one-step", the function will skip resampling the datasets, and directly resample from the original list
 ### of ExpressionSets 
 balance.variables=NULL
 ### balance.variables will be a vector of covariate names that should be balanced in the simulation. After balancing, the prevalence of covariate in each 
 ### result set should be the same as the overall distribution across all original data sets. Default is set as NULL, when it will not balance over any covariate.
){  
  #### Step.1 Drawing sets
  if(type !="one-step" && type != "two-steps")
    stop("Wrong type.")
  else if(type == "one-step")
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
  new.y.vars <- list()
  if(!is.null(balance.variables)){      
    print(paste("covariate: ", balance.variables, sep=""))
    for(i in 1:length(esets)){
      sampleind[[i]] <- sample(1:length(sampleNames(esets[[setsID[i]]])), n.samples, replace=TRUE, prob=probs.list[[setsID[i]]])
      samplesets[[i]] <- esets[[setsID[i]]][, sampleind[[i]]]
      if(length(y.vars)!=0){
        new.y.vars[[i]] <- y.vars[[setsID[i]]][sampleind[[i]], ]
      } 
    }        
  }    
  else {
    print("covariate: NULL")
    for(i in 1:length(esets)){
      sampleind[[i]] <- sample(1:length(sampleNames(esets[[setsID[i]]])), n.samples, replace=TRUE)
      samplesets[[i]] <- esets[[setsID[i]]][, sampleind[[i]]]
      if(length(y.vars)!=0){
        new.y.vars[[i]] <- y.vars[[setsID[i]]][sampleind[[i]], ]
      }
    }
  }
  if(length(y.vars)==0) new.y.vars <- y.vars
  names(samplesets) <- setsID
  ## We still need to remember which studies these samples are from, because we will need to 
  ## generate survival time later, and beta from CoxBoost should be relevant to original data sets 
  ## End Step.2
  res <- list(esets=samplesets, indices=sampleind, setsID=setsID, y.vars=new.y.vars, 
              prob.desired=prob_desired, prob.real=probs.list)
  return(res)
  ### returns a list of simulated ExpressionSets, with names indicating its original set, and indices of the original patients.
},ex=function(){
  library(curatedOvarianData)
  data( E.MTAB.386_eset )
  eset1 <- E.MTAB.386_eset[1:10, 1:5]
  eset2 <- E.MTAB.386_eset[1:10, 6:10]
  eset3 <- E.MTAB.386_eset[1:10, 11:15]
    
  ## simulate on multiple ExpressionSets
  esets.list <- list(eset1, eset2, eset3)  
  # one-step bootstrap: skip resampling set labels
  simmodels <- simData(esets.list, 20, type="one-step")  
  # two-step-non-parametric bootstrap
  simmodels <- simData(esets.list, 10, type="two-steps")
  
  ## simulate one set
  simmodels <- simData(list(eset1), 10, type="two-steps")
  
  ## balancing covariates
  # single covariate
  simmodels <- simData(list(eset2), 5, balance.variables="tumorstage")
  # check the balancing effect
  simmodels$prob.desired
  table(eset2$tumorstage) / sum(table(eset2$tumorstage))
  table(simmodels$esets[[1]]$tumorstage) / sum( table(simmodels$esets[[1]]$tumorstage))
  # multiple covariates
  simmodels <- simData(list(eset2), 5, 
                       balance.variables=c("tumorstage", "age_at_initial_pathologic_diagnosis"))  
})