funCV <- function
### Cross validation function
(eset,
 ### simulated ExpressionSets 
 fold, 
 ### tells the number of folds in cross validation
 p.cutoff, 
 ### trainFun parameter, in this case passes to masomenos()
 trainFun,
 ### training function, in this case masomenos()
 testFun, 
 ### test function, in this case funTest()
 cstatFun
 ### function for calculating the c stat, in this case calcu_cstat()
){
  print("CV")
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
  ## cross validation
  z <- numeric(fold)
  for(i in 1:fold){
    print(paste("fold = ", i, sep=""))
    testeset <- eset[, setindex[[i]]]
    traineset <- eset[, -(setindex[[i]])]
    trainX <- t(exprs(traineset))
    trainY <- pData(traineset)[, c("dmfs.time","dmfs.cens")]
    trainY[, 1] <- as.numeric(as.character(trainY[, 1]))
    trainY[, 2] <- as.numeric(as.character(trainY[, 2]))
    beta <- trainFun(trainX, trainY, length(featureNames(eset)), p.cutoff)
    z[i] <- testFun(eset=testeset, masomenosBeta=beta, cstatFun=cstatFun)
  }
  z <- sum(z) / fold
  return(z)
  ### returns the c statistics of cross validation(CV)
}