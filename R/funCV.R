funCV <- structure(function
### Cross validation function
(eset,
 ### one ExpressionSet to do cross validation with 
 fold,  
 ### the number of folds in cross validation
 y.var,
 ### response variable, matrix, data.frame(with 2 columns) or Surv object
 trainFun=masomenos,
 ### training function, which takes gene expression matrix X and response variable y as input, the coefficients as output 
 funCvSubset=cvSubsets
 ### function to divide one Expression Set into subsets for cross validation
){
  setindex <- funCvSubset(eset, fold) 
  ## cross validation
  z <- numeric(fold)
  for(i in 1:fold){
    print(paste("fold = ", i, sep=""))
    testeset <- eset[, setindex[[i]]]
    testX <- t(exprs(testeset))
    testY <- y.var[setindex[[i]], ]
    testY[, 1] <- as.numeric(as.character(testY[, 1]))
    testY[, 2] <- as.numeric(as.character(testY[, 2]))
    testY <- Surv(testY[, 1], testY[, 2])
    
    traineset <- eset[, -(setindex[[i]])]
    trainX <- t(exprs(traineset))
    trainY <- y.var[-setindex[[i]], ]
    trainY[, 1] <- as.numeric(as.character(trainY[, 1]))
    trainY[, 2] <- as.numeric(as.character(trainY[, 2]))
    trainY <- Surv(trainY[, 1], trainY[, 2])
    
    beta <- trainFun(trainX, trainY)
    lp <- testX %*% beta    
	  z[i] <- rcorr.cens(-lp, testY)[1]
  }
  z <- sum(z) / fold
  return(z)
  ### returns the c statistics of cross validation(CV)
},ex=function(){
  library(curatedOvarianData)
  data( E.MTAB.386_eset )
  eset <- E.MTAB.386_eset[1:100, 1:30]
  time <- sample(2000:5000, 30)
  cens <- sample(0:1, 30, replace=TRUE)
  y <- Surv(time, cens)  
  y1 <- cbind(time, cens)
  
  z1.score <- funCV(eset, 3, y)
  z2.score <- funCV(eset, 3, y1, trainFun=plusMinus)
  ## any training function will do as long as it takes the gene expression matrix X
  ## and response variable y(matrix, data.frame or Surv object) as parameters, and
  ## return the coefficients as its value
})