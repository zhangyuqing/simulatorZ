simData <- structure(function
### simData is a function to perform non-parametric bootstrap resampling
### on a list of (original) data sets, both on set level and patient level,
### in order to simulate independent genomic sets. 
(obj,
 ### a list of ExpressionSets, matrics or SummarizedExperiments. If elements are 
 ### matricse, columns represent samples
 n.samples,
 ### an integer indicating how many samples should be resampled from each set
 y.vars=list(),
 ### a list of response variables, can be Surv object, or matrix or data.frame
 ### with two columns
 type="two-steps",
 ### string "one-step" or "two-steps". If type="one-step", the function will 
 ### skip resampling the datasets, and directly resample from the original list
 ### of obj 
 balance.variables=NULL
 ### balance.variables will be a vector of covariate names that should be 
 ### balanced in the simulation. After balancing, the prevalence of covariate 
 ### in each result set should be the same as the overall distribution across 
 ### all original data sets. Default is set as NULL, when it will not balance 
 ### over any covariate. if isn't NULL, esets parameter should only be of class
 ### ExpressionSet
 ){  
  #### Step.1 Drawing sets
  if(type !="one-step" && type != "two-steps"){
    stop("Wrong type.")
  }    
  else if(type == "one-step"){
    setsID <- 1:length(obj)
  }    
  else{
    prob.set <- rep((1/length(obj)), times=length(obj))
    setsID <- sample(1:length(obj), prob=prob.set, replace=TRUE)  # labels of data sets
  }  
  print(setsID) 
  
  if(class(obj[[1]])=="ExpressionSet"){
    #### Step.2 Calculate the joint probability distribution  
    ## rbinding all covariates data  
    covariates.list <- lapply(obj, function(obj.ele){
      if(length(balance.variables) == 1){
        return(as.character(pData(obj.ele)[, balance.variables]))
      }
      else{
        return(as.character(do.call(paste, pData(obj.ele)[, balance.variables])))
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
  } 
  
  #### Step. 3 Drawing patients
  samplesets <- sampleind <- list()
  new.y.vars <- list()
  
  
  
  if(!is.null(balance.variables)){
    if(class(obj[[1]])!="ExpressionSet"){
      stop("Cannot balance covariates for type other than ExpressionSet!")
    }
    else{
      print(paste("covariate: ", balance.variables, sep=""))
      for(i in 1:length(obj)){        
        sampleind[[i]] <- sample(1:length(sampleNames(obj[[setsID[i]]])), n.samples, replace=TRUE, prob=probs.list[[setsID[i]]])
        samplesets[[i]] <- obj[[setsID[i]]][, sampleind[[i]]]
        if(length(y.vars)!=0){
          new.y.vars[[i]] <- y.vars[[setsID[i]]][sampleind[[i]], ]
        } 
      }
    }            
  }    
  else {
    print("covariate: NULL")
    for(i in 1:length(obj)){
      if(class(obj[[1]])=="ExpressionSet")
        num.sam <- ncol(exprs(obj[[setsID[i]]]))
      else if(class(obj[[1]])=="matrix")
        num.sam <- ncol(obj[[setsID[i]]])
      else if(class(obj[[1]])=="SummarizedExperiment")
        num.sam <- ncol(assay(obj[[setsID[i]]]))
      sampleind[[i]] <- sample(1:num.sam, n.samples, replace=TRUE)
      samplesets[[i]] <- obj[[setsID[i]]][, sampleind[[i]]]
      if(length(y.vars)!=0){
        new.y.vars[[i]] <- y.vars[[setsID[i]]][sampleind[[i]], ]
      }
    }
    prob_desired <- c()
    probs.list <- list()
  }
  if(length(y.vars)==0) new.y.vars <- y.vars
  names(samplesets) <- setsID
  ## We still need to remember which studies these samples are from, because we will need to 
  ## generate survival time later, and beta from CoxBoost should be relevant to original data sets 
  ## End Step.2
  res <- list(obj=samplesets, indices=sampleind, setsID=setsID, y.vars=new.y.vars, 
              prob.desired=prob_desired, prob.real=probs.list)
  return(res)
  ### returns a list of simulated ExpressionSets, with names indicating its original set, and indices of the original patients.
  ### prob.desired and prob.real are only useful when balance.varaibles is set.
  ### prob.desired shows overall distrubition of the specified covariate. prob.list
  ### shows the sampling probability in each set after balancing
},ex=function(){
  library(curatedOvarianData)
  library(GenomicRanges)
  source(system.file("extdata", "patientselection.config",
                     package="curatedOvarianData"))
  source(system.file("extdata", "createEsetList.R", package="curatedOvarianData"))
  esets.list <- lapply(esets, function(eset){
    return(eset[1:1000, 1:10])
  })
  
  ## simulate on multiple ExpressionSets
  set.seed(8)
  # one-step bootstrap: skip resampling set labels
  simmodels <- simData(esets.list, 20, type="one-step")  
  # two-step-non-parametric bootstrap
  simmodels <- simData(esets.list, 10, type="two-steps")
  
  ## simulate one set
  simmodels <- simData(list(esets.list[[1]]), 10, type="two-steps")
  
  ## balancing covariates
  # single covariate
  simmodels <- simData(list(esets.list[[1]]), 5, balance.variables="tumorstage")
  
  # multiple covariates
  simmodels <- simData(list(esets.list[[1]]), 5, 
                       balance.variables=c("tumorstage", "age_at_initial_pathologic_diagnosis"))  

  ## Support matrices
  X.list <- lapply(esets.list, function(eset){
    return(exprs(eset))
  })
  simmodels <- simData(X.list, 20, type="two-steps")
  
  ## Support SummarizedExperiment
  nrows <- 200; ncols <- 6
  counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
  rowData <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),
                     IRanges(floor(runif(200, 1e5, 1e6)), width=100),
                     strand=sample(c("+", "-"), 200, TRUE))
  colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),
                       row.names=LETTERS[1:6])
  sset <- SummarizedExperiment(assays=SimpleList(counts=counts),
                               rowData=rowData, colData=colData)
  
  s.list <- list(sset[,1:3], sset[,4:6])
  simmodels <- simData(s.list, 20, type="two-steps")
})