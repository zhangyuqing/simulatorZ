zmatrix <- function
### Takes a list of ExpressionSets to generate one matrix of validation statistics
### already-defined training, test and CV functions. 
(esets, 
 ### simulated ExpressionSets 
 y.vars, 
 ### some strings to indicate the response variable
 trainingFun, 
 ### training function, in this case masomenos()
 testFun, 
 ### test function, in this case funTest()
 cvFun, 
 ### function to do cross validation, in this case funCV()
 cstatFun,
 ### function to calculate the c statistics, in this case calcu_cstat
 p.cutoff, 
 ### trainFun parameter, in this case passes to masomenos()
 fold
 ### cvFun parameter, in this case passes to funCV()
){
  Zmatrix <- matrix(, nrow=length(esets), ncol=length(esets))
  for(i in 1:length(esets)){
    for(j in 1:length(esets)){
      if(i == j) {
        print(i)
        Zmatrix[i, j] <- cvFun(eset=(esets[[i]]), fold, p.cutoff, 
                               trainingFun, testFun, cstatFun)
        print(Zmatrix[i, j])
      }
      else if (i != j){
        print(paste("train:", i, "test:", j, sep=" "))
        trainX <- t(exprs((esets[[i]])))
        trainX <- scale(trainX)
        trainY <- pData(esets[[i]])[, y.vars]
        trainY[, 1] <- as.numeric(as.character(trainY[, 1]))
        trainY[, 2] <- as.numeric(as.character(trainY[, 2]))
        beta <- trainingFun(trainX, trainY, length(featureNames(esets[[i]])), p.cutoff)
        Zmatrix[i, j] <- testFun(esets[[j]], beta, cstatFun)
        print(Zmatrix[i, j])
      }
    }
  }
  return(Zmatrix)
  ### outputs one matrix of validation statistics
}