zmatrix <- structure(function
### generate a matrix of c statistics
(esets, 
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
  Zmatrix <- matrix(, nrow=length(esets), ncol=length(esets))
  for(i in 1:length(esets)){
    for(j in 1:length(esets)){
      if(i == j) {
        print(paste("CV: ", i, sep=""))
        Zmatrix[i, j] <- cvFun(eset=esets[[i]], fold=fold, y.var=y.vars[[i]],
							   trainFun=trainingFun, funCvSubset=cvSubsetFun)
        print(Zmatrix[i, j])
      }
      else if (i != j){
        print(paste("train:", i, "test:", j, sep=" "))
        
        if(class(esets[[1]])=="ExpressionSet"){
          testX <- t(exprs(esets[[j]]))
          trainX <- t(exprs(esets[[i]]))
        }
        else if(class(esets[[1]])=="matrix"){
          testX <- t(esets[[j]])
          trainX <- t(esets[[i]])
        }
        else if(class(esets[[1]])=="SummarizedExperiment"){
          testX <- t(assay(esets[[j]]))
          trainX <- t(assay(esets[[i]]))
        }
        else stop("Wrong class of eset!")
        
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
  data( E.MTAB.386_eset )
  eset1 <- E.MTAB.386_eset[1:100, 1:30]
  eset2 <- E.MTAB.386_eset[1:100, 31:60]
  eset3 <- E.MTAB.386_eset[1:100, 61:90]  
  esets <- list(eset1, eset2, eset3) 
  
  time1 <- eset1$days_to_death
  #cens1 <- c(0, 0, 0, 1, 1)
  cens1 <- sample(0:1, 30, replace=TRUE)
  y1 <- Surv(time1, cens1)
  time2 <- eset2$days_to_death
  #cens2 <- c(1, 1, 0, 0, 0)
  cens2 <- sample(0:1, 30, replace=TRUE)
  y2 <- Surv(time2, cens2)
  time3 <- eset3$days_to_death
  #cens3 <- c(1, 0, 0, 0, 1)
  cens3 <- sample(0:1, 30, replace=TRUE)
  y3 <- Surv(time3, cens3)
  y.vars <- list(y1, y2, y3)
  
  # generate on original ExpressionSets
  z <- zmatrix(esets, y.vars, 3)
  
  # generate on simulated ExpressionSets
  simmodels <- simBootstrap(esets, y.vars, 10, 100)
  z <- zmatrix(simmodels$esets.list, simmodels$y.vars.list, 3)
  
  # support matrix
  X.list <- lapply(esets, function(eset){
    newx <- exprs(eset) ### columns represent samples !!
    return(newx)
  }) 
  z <- zmatrix(X.list, y.vars, 3)
  
  # support SummarizedExperiment
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
  
  z <- zmatrix(list(sset[,1:3], sset[,4:6]), list(y.vars[1:3,],y.vars[4:6,]), 3)
})