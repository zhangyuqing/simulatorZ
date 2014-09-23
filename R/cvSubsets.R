cvSubsets <- structure(function
### To generate a list of subsets(indices of observations) from one set
(obj,
### a ExpressionSet, matrix or SummarizedExperiment object. If it is a matrix,
### columns represent samples
fold
### the number of folds in cross validation. 
### Number of observations in the set does not need to be a multiple of fold
){
  setindex <- list()
  if(class(obj)=="ExpressionSet")
    seq <- 1:ncol(exprs(obj))
  else if(class(obj)=="matrix")
    seq <- 1:ncol(obj)
  else if(class(obj)=="SummarizedExperiment")
    seq <- 1:ncol(assay(obj))
  else stop("Wrong class of obj!")
    
  n.subset <- round(length(seq) / fold)
  ## divide into 4 subsets randomly   
  for(i in 1:(fold-1)){
    if(class(obj)=="ExpressionSet")
      seq <- 1:ncol(exprs(obj))
    else if(class(obj)=="matrix")
      seq <- 1:ncol(obj)
    else if(class(obj)=="SummarizedExperiment")
      seq <- 1:ncol(assay(obj))
    if(length(do.call(c, setindex)) > 0) seq <- seq[-do.call(c, setindex)]
    setindex[[i]] <- sample(seq, n.subset, replace=FALSE)  
  }
  if(class(obj)=="ExpressionSet")
    seq <- 1:ncol(exprs(obj))
  else if(class(obj)=="matrix")
    seq <- 1:ncol(obj)
  else if(class(obj)=="SummarizedExperiment")
    seq <- 1:ncol(assay(obj))
  setindex[[fold]] <- seq[-do.call(c, setindex)]
  return(setindex)
  ### returns the list of indices of subsets
},ex=function(){
  library(curatedOvarianData)
  data(E.MTAB.386_eset)
  
  set.seed(8)
  id <- cvSubsets(E.MTAB.386_eset, 3)
  
  subset1 <- E.MTAB.386_eset[, id[[1]]]
  subset2 <- E.MTAB.386_eset[, id[[2]]]
  subset3 <- E.MTAB.386_eset[, id[[3]]]
  
  ## Number of observations in the set does not need to be a multiple of
  ## the fold parameter
  id2 <- cvSubsets(E.MTAB.386_eset, 5)
  subsets <- list()
  subsets[[1]] <- E.MTAB.386_eset[, id2[[1]]]
  subsets[[2]] <- E.MTAB.386_eset[, id2[[2]]]
  subsets[[3]] <- E.MTAB.386_eset[, id2[[3]]]
  subsets[[4]] <- E.MTAB.386_eset[, id2[[4]]]
  subsets[[5]] <- E.MTAB.386_eset[, id2[[5]]]
})