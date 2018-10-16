cvSubsets <- structure(function
### To generate a list of subsets(indices of observations) from one set
(obj,
### a ExpressionSet, matrix or RangedSummarizedExperiment object. If it is a
### matrix,columns represent samples
fold
### the number of folds in cross validation. 
### Number of observations in the set does not need to be a multiple of fold
){
  setindex <- list()
  obj <- getMatrix(obj)
  n.subset <- round(ncol(obj) / fold)
  ## divide into 4 subsets randomly   
  for(i in 1:(fold-1)){
    seq <- seq_len(ncol(obj))
    if(length(do.call(c, setindex)) > 0) seq <- seq[-do.call(c, setindex)]
    setindex[[i]] <- sample(seq, n.subset, replace=FALSE)  
  }
  seq <- seq_len(ncol(obj))
  setindex[[fold]] <- seq[-do.call(c, setindex)]
  return(setindex)
  ### returns the list of indices of subsets
},ex=function(){
  library(curatedOvarianData)
  data(E.MTAB.386_eset)
  
  id <- cvSubsets(E.MTAB.386_eset, 3)
  subsets <- lapply(1:3, function(i){E.MTAB.386_eset[1:10, id[[i]]]})
  sapply(subsets, dim)
  rm(subsets)
  
  ## Number of observations in the set does not need to be a multiple of
  ## the fold parameter
  id2 <- cvSubsets(E.MTAB.386_eset, 5)
  subsets <- lapply(1:5, function(j){E.MTAB.386_eset[1:10, id2[[j]]]})
  sapply(subsets, dim)
  rm(subsets)
})
