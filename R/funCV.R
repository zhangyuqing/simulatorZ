funCV <- structure(function
### Cross validation function
(obj,
 ### a ExpressionSet, matrix or RangedSummarizedExperiment object. If it is a
 ### matrix, columns represent samples 
 fold,  
 ### the number of folds in cross validation
 y.var,
 ### response variable, matrix, data.frame(with 2 columns) or Surv object
 trainFun=masomenos,
 ### training function, which takes gene expression matrix X and response variable y as input, the coefficients as output 
 funCvSubset=cvSubsets,
 ### function to divide one Expression Set into subsets for cross validation
 covar=NULL
 ### other covariates to be added in as predictors
){
  setindex <- funCvSubset(obj, fold)

  if(!is.null(covar)){
    X1 <- pData(obj)[,covar]
    X2 <- changevar(X1, covar)
    obj <- t(getMatrix(obj))
    obj <- t(cbind(X2, obj))
  }else{
    obj <- getMatrix(obj)
  }
  
  ## cross validation
  z <- numeric(fold)
  for(i in 1:fold){
    print(paste("fold = ", i, sep=""))
    testX <- t(obj[, setindex[[i]]])
    trainX <- t(obj[, -(setindex[[i]])])
    testY <- y.var[setindex[[i]], ]
    testY[, 1] <- as.numeric(as.character(testY[, 1]))
    testY[, 2] <- as.numeric(as.character(testY[, 2]))
    testY <- Surv(testY[, 1], testY[, 2])    
    trainY <- y.var[-setindex[[i]], ]
    trainY[, 1] <- as.numeric(as.character(trainY[, 1]))
    trainY[, 2] <- as.numeric(as.character(trainY[, 2]))
    trainY <- Surv(trainY[, 1], trainY[, 2])
    
    # handle duplicate sample names if there is bootstrap resampling
    if(all(!duplicated(rownames(trainX)))){rownames(trainX) <- NULL}
    
    beta <- trainFun(trainX, trainY)
    lp <- testX %*% beta    
    z[i] <- rcorr.cens(-lp, testY)["C Index"]
  }
  z <- sum(z,na.rm=TRUE) / sum(!is.nan(z))
  return(z)
  ### returns the c statistics of cross validation(CV)
},ex=function(){
  library(curatedOvarianData)
  library(GenomicRanges)
  set.seed(8)
  data( E.MTAB.386_eset )
  eset <- E.MTAB.386_eset[1:100, 1:30]
  rm(E.MTAB.386_eset)
  
  time <- eset$days_to_death
  cens.chr <- eset$vital_status
  cens <- rep(0, length(cens.chr))
  cens[cens.chr=="living"] <- 1
  y <- Surv(time, cens)  
  y1 <- cbind(time, cens)
  
  nrows <- 200; ncols <- 6
  counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
  rowRanges <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),
                       IRanges(floor(runif(200, 1e5, 1e6)), width=100),
                       strand=sample(c("+", "-"), 200, TRUE))
  colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),
                       row.names=LETTERS[1:6])
  sset <- SummarizedExperiment(assays=SimpleList(counts=counts),
                               rowRanges=rowRanges, colData=colData)
  time <- c(1588,1929,1813,1542,1830,1775)  
  cens <- c(1,0,1,1,1,1)
  y.vars <- Surv(time, cens)
  
  funCV(eset, 3, y)
  funCV(exprs(eset), 3, y1)
  funCV(sset, 3, y.vars)
  ## any training function will do as long as it takes the gene expression matrix X
  ## and response variable y(matrix, data.frame or Surv object) as parameters, and
  ## return the coefficients as its value
})
