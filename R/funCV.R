funCV <- structure(function
### Cross validation function
(obj,
 ### a ExpressionSet, matrix or SummarizedExperiment object. If it is a matrix,
 ### columns represent samples 
 fold,  
 ### the number of folds in cross validation
 y.var,
 ### response variable, matrix, data.frame(with 2 columns) or Surv object
 trainFun=masomenos,
 ### training function, which takes gene expression matrix X and response variable y as input, the coefficients as output 
 funCvSubset=cvSubsets
 ### function to divide one Expression Set into subsets for cross validation
){
  setindex <- funCvSubset(obj, fold) 
  ## cross validation
  z <- numeric(fold)
  for(i in 1:fold){
    print(paste("fold = ", i, sep=""))
    testeset <- obj[, setindex[[i]]]
    traineset <- obj[, -(setindex[[i]])]
    
    if(class(obj)=="ExpressionSet"){
      testX <- t(exprs(testeset))
      trainX <- t(exprs(traineset))
    }
    else if(class(obj)=="matrix"){
      testX <- t(testeset)
      trainX <- t(traineset)
    }
    else if(class(obj)=="SummarizedExperiment"){
      testX <- t(assay(testeset))
      trainX <- t(assay(traineset))
    }
    else stop("Wrong class of obj!")
          
    testY <- y.var[setindex[[i]], ]
    testY[, 1] <- as.numeric(as.character(testY[, 1]))
    testY[, 2] <- as.numeric(as.character(testY[, 2]))
    testY <- Surv(testY[, 1], testY[, 2])    
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
  library(GenomicRanges)
  set.seed(8)
  data( E.MTAB.386_eset )
  eset <- E.MTAB.386_eset[1:100, 1:30]
  time <- eset$days_to_death
  cens.chr <- eset$vital_status
  cens <- c()
  for(i in 1:length(cens.chr)){
    if(cens.chr=="living") cens[i] <- 1
    else cens[i] <- 0
  }
  y <- Surv(time, cens)  
  y1 <- cbind(time, cens)
  
  nrows <- 200; ncols <- 6
  counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
  rowData <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),
                     IRanges(floor(runif(200, 1e5, 1e6)), width=100),
                     strand=sample(c("+", "-"), 200, TRUE))
  colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),
                       row.names=LETTERS[1:6])
  sset <- SummarizedExperiment(assays=SimpleList(counts=counts),
                               rowData=rowData, colData=colData)
  time <- sample(4500:4700, 6, replace=TRUE)
  cens <- sample(0:1, 6, replace=TRUE)
  y.vars <- Surv(time, cens)
  
  funCV(eset, 3, y)
  funCV(eset, 3, y1, trainFun=plusMinus)
  funCV(exprs(eset), 3, y)
  
  funCV(sset, 3, y.vars)
  ## any training function will do as long as it takes the gene expression matrix X
  ## and response variable y(matrix, data.frame or Surv object) as parameters, and
  ## return the coefficients as its value
})