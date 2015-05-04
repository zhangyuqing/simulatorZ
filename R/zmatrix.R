zmatrix <- structure(function
### generate a matrix of c statistics
(obj, 
 ### a list of ExpressionSet, matrix or SummarizedExperiment objects, if its 
 ### elements are matrices, columns represent samples 
 y.vars, 
 ### a list of response variables, all the response variables shold be
 ### matrix, data.frame(with 2 columns) or Surv object
 fold,
 ### cvFun parameter, in this case passes to funCV()
 trainingFun=masomenos, 
 ### training function
 cvFun=funCV, 
 ### function to perform cross study within one set
 cvSubsetFun=cvSubsets
 ### function to divide the expression sets into subsets for cross validation 
){
  Zmatrix <- matrix(, nrow=length(obj), ncol=length(obj))
  for(i in 1:length(obj)){
    for(j in 1:length(obj)){
      if(i == j) {
        print(paste("CV: ", i, sep=""))
        Zmatrix[i, j] <- cvFun(obj=obj[[i]], fold=fold, y.var=y.vars[[i]],
							   trainFun=trainingFun, funCvSubset=cvSubsetFun)
        print(Zmatrix[i, j])
      }
      else if (i != j){
        print(paste("train:", i, "test:", j, sep=" "))
        
        if(class(obj[[1]])=="ExpressionSet"){
          testX <- t(exprs(obj[[j]]))
          trainX <- t(exprs(obj[[i]]))
        }
        else if(class(obj[[1]])=="matrix"){
          testX <- t(obj[[j]])
          trainX <- t(obj[[i]])
        }
        else if(class(obj[[1]])=="SummarizedExperiment"){
          testX <- t(assay(obj[[j]]))
          trainX <- t(assay(obj[[i]]))
        }
        else stop("Wrong class of object!")
        
        trainY <- y.vars[[i]]
        trainY[, 1] <- as.numeric(as.character(trainY[, 1]))
        trainY[, 2] <- as.numeric(as.character(trainY[, 2]))
        trainY <- Surv(trainY[, 1], trainY[, 2])
	    	testY <- y.vars[[j]]
        testY[, 1] <- as.numeric(as.character(testY[, 1]))
        testY[, 2] <- as.numeric(as.character(testY[, 2]))
	    	testY <- Surv(testY[, 1], testY[, 2])
        
        beta <- trainingFun(trainX, trainY)
		    lp <- testX %*% beta
        Zmatrix[i, j] <- rcorr.cens(-lp, testY)[1]
        print(Zmatrix[i, j])
      }
    }
  }
  return(Zmatrix)
  ### outputs one matrix of validation statistics
},ex=function(){
  library(curatedOvarianData)
  library(GenomicRanges)
  data(GSE17260_eset)
  data(E.MTAB.386_eset)
  data(GSE14764_eset)
  esets <- list(GSE17260=GSE17260_eset, E.MTAB.386=E.MTAB.386_eset, GSE14764=GSE14764_eset)
  esets.list <- lapply(esets, function(eset){
    return(eset[1:500, 1:20])
  })
  
  ## simulate on multiple ExpressionSets
  set.seed(8) 
  
  y.list <- lapply(esets.list, function(eset){
    time <- eset$days_to_death
    cens.chr <- eset$vital_status
    cens <- c()
    for(i in 1:length(cens.chr)){
      if(cens.chr[i] == "living") cens[i] <- 1
      else cens[i] <- 0
    }
    y <- Surv(time, cens)
    return(y)
  })
  
  # generate on original ExpressionSets
  z <- zmatrix(esets.list, y.list, 3)
  
  # generate on simulated ExpressionSets
  simmodels <- simBootstrap(esets.list, y.list, 10, 100)
  z <- zmatrix(simmodels$obj.list, simmodels$y.vars.list, 3)
  
  # support matrix
  X.list <- lapply(esets.list, function(eset){
    newx <- exprs(eset) ### columns represent samples !!
    return(newx)
  }) 
  z <- zmatrix(X.list, y.list, 3)
  
  # support SummarizedExperiment
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
