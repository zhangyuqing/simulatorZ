zmatrix <- structure(function
### generate a matrix of c statistics
(obj, 
 ### a list of ExpressionSet, matrix or RangedSummarizedExperiment objects.
 ### If its elements are matrices, columns represent samples 
 y.vars, 
 ### a list of response variables, all the response variables shold be
 ### matrix, data.frame(with 2 columns) or Surv object
 fold,
 ### cvFun parameter, in this case passes to funCV()
 trainingFun=masomenos, 
 ### training function
 cvFun=funCV, 
 ### function to perform cross study within one set
 cvSubsetFun=cvSubsets,
 ### function to divide the expression sets into subsets for cross validation 
 covar=NULL
 ### other covariates to be added as predictors
){
  Zmatrix <- matrix(, nrow=length(obj), ncol=length(obj))
  for(i in seq_along(obj)){
    for(j in seq_along(obj)){
      if(i == j) {
        print(paste("CV: ", i, sep=""))
        Zmatrix[i, j] <- cvFun(obj=obj[[i]], fold=fold, y.var=y.vars[[i]],
                               trainFun=trainingFun, funCvSubset=cvSubsetFun,
                               covar=covar)
        print(Zmatrix[i, j])
      }
      else if (i != j){
        print(paste("train:", i, "test:", j, sep=" "))
        trainX <- t(getMatrix(obj[[i]]))
        testX <- t(getMatrix(obj[[j]]))
        
        if(!is.null(covar)){
          X1 <- pData(obj[[i]])[, covar]
          X2 <- changevar(X1, covar)
          trainX <- cbind(X2, trainX)
          
          X3 <- pData(obj[[j]])[, covar]
          X4 <- changevar(X3, covar)
          testX <- cbind(X4, testX)
        }
        
        trainY <- y.vars[[i]]
        trainY[, 1] <- as.numeric(as.character(trainY[, 1]))
        trainY[, 2] <- as.numeric(as.character(trainY[, 2]))
        trainY <- Surv(trainY[, 1], trainY[, 2])
        testY <- y.vars[[j]]
        testY[, 1] <- as.numeric(as.character(testY[, 1]))
        testY[, 2] <- as.numeric(as.character(testY[, 2]))
        testY <- Surv(testY[, 1], testY[, 2])
        
        # handle duplicate sample names if there is bootstrap resampling
        if(all(!duplicated(rownames(trainX)))){rownames(trainX) <- NULL}
        
        beta <- trainingFun(trainX, trainY)
        lp <- testX %*% beta
        Zmatrix[i, j] <- rcorr.cens(-lp, testY)["C Index"]
        print(Zmatrix[i, j])
      }
    }
  }
  return(Zmatrix)
  ### outputs one matrix of validation statistics
},ex=function(){
  library(curatedOvarianData)
  library(GenomicRanges)
  data(E.MTAB.386_eset)
  data(GSE14764_eset)
  esets.list <- list(E.MTAB.386=E.MTAB.386_eset[1:100, 1:30], GSE14764=GSE14764_eset[1:100, 1:30])
  rm(E.MTAB.386_eset, GSE14764_eset)
  
  ## simulate on multiple ExpressionSets
  set.seed(8) 
  
  y.list <- lapply(esets.list, function(eset){
    time <- eset$days_to_death
    cens.chr <- eset$vital_status
    cens <- rep(0, length(cens.chr))
    cens[cens.chr=="living"] <- 1
    return(Surv(time, cens))
  })
  
  # generate on original ExpressionSets
  z <- zmatrix(esets.list, y.list, 3)
  
  # generate on simulated ExpressionSets
  simmodels <- simBootstrap(esets.list, y.list, 100, 100)
  z <- zmatrix(simmodels$obj.list, simmodels$y.vars.list, 3)
  
  # support matrix
  X.list <- lapply(esets.list, function(eset){
    return(exprs(eset)) ### columns represent samples !!
  }) 
  z <- zmatrix(X.list, y.list, 3)
  
  # support RangedSummarizedExperiment
  nrows <- 200; ncols <- 6
  counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
  rowRanges <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),
                       IRanges(floor(runif(200, 1e5, 1e6)), width=100),
                       strand=sample(c("+", "-"), 200, TRUE))
  colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),
                       row.names=LETTERS[1:6])
  sset <- SummarizedExperiment(assays=SimpleList(counts=counts),
                               rowRanges=rowRanges, colData=colData)
  
  time <- sample(4500:4700, 6, replace=TRUE)
  cens <- sample(0:1, 6, replace=TRUE)
  y.vars <- Surv(time, cens)
  
  z <- zmatrix(list(sset[,1:3], sset[,4:6]), list(y.vars[1:3,],y.vars[4:6,]), 3)
})
