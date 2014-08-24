cvSubsets <- structure(function
### To generate a list of subsets(indices of observations) from one Expression Set, to do cross validation
(eset,
### one ExpressionSet to do cross validation with 
fold
### the number of folds in cross validation
){
  setindex <- list()
  seq <- 1:length(sampleNames(eset))
  n.subset <- round(length(seq) / fold)
  ## divide into 4 subsets randomly   
  for(i in 1:(fold-1)){
    seq <- 1:length(sampleNames(eset))
    if(length(do.call(c, setindex)) > 0) seq <- seq[-do.call(c, setindex)]
    setindex[[i]] <- sample(seq, n.subset, replace=FALSE)  
  }
  seq <- 1:length(sampleNames(eset))
  setindex[[fold]] <- seq[-do.call(c, setindex)]
  return(setindex)
  ### returns the list of indices of subsets
},ex=function(){
  library(curatedOvarianData)
  data(E.MTAB.386_eset)
  id <- cvSubsets(E.MTAB.386_eset, 3)
  
  subset1 <- E.MTAB.386_eset[, id[[1]]]
  subset2 <- E.MTAB.386_eset[, id[[2]]]
  subset3 <- E.MTAB.386_eset[, id[[3]]]
  subset1
  
  ## Number of observations in the set does not need to be a multiple of
  ## the fold parameter
  id2 <- cvSubsets(E.MTAB.386_eset, 5)
  subsets <- list()
  subsets[[1]] <- E.MTAB.386_eset[, id2[[1]]]
  subsets[[2]] <- E.MTAB.386_eset[, id2[[2]]]
  subsets[[3]] <- E.MTAB.386_eset[, id2[[3]]]
  subsets[[4]] <- E.MTAB.386_eset[, id2[[4]]]
  subsets[[5]] <- E.MTAB.386_eset[, id2[[5]]]
  subsets
})