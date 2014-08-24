cvSubsets <- structure(function
### To generate a list of subsets(indices of observations) from one set
(eset,
### a ExpressionSet, matrix or SummarizedExperiment object. If it is a matrix,
### columns represent samples
fold
### the number of folds in cross validation. 
### Number of observations in the set does not need to be a multiple of fold
){
  setindex <- list()
  if(class(eset)=="ExpressionSet")
    seq <- 1:ncol(exprs(eset))
  else if(class(eset)=="matrix")
    seq <- 1:ncol(eset)
  else if(class(eset)=="SummarizedExperiment")
    seq <- 1:ncol(assay(eset))
  else stop("Wrong class of eset!")
    
  n.subset <- round(length(seq) / fold)
  ## divide into 4 subsets randomly   
  for(i in 1:(fold-1)){
    if(class(eset)=="ExpressionSet")
      seq <- 1:ncol(exprs(eset))
    else if(class(eset)=="matrix")
      seq <- 1:ncol(eset)
    else if(class(eset)=="SummarizedExperiment")
      seq <- 1:ncol(assay(eset))
    if(length(do.call(c, setindex)) > 0) seq <- seq[-do.call(c, setindex)]
    setindex[[i]] <- sample(seq, n.subset, replace=FALSE)  
  }
  if(class(eset)=="ExpressionSet")
    seq <- 1:ncol(exprs(eset))
  else if(class(eset)=="matrix")
    seq <- 1:ncol(eset)
  else if(class(eset)=="SummarizedExperiment")
    seq <- 1:ncol(assay(eset))
  setindex[[fold]] <- seq[-do.call(c, setindex)]
  return(setindex)
  ### returns the list of indices of subsets
},ex=function(){
  library(curatedOvarianData)
  library(GenomicRanges)
  data(E.MTAB.386_eset)
  id <- cvSubsets(E.MTAB.386_eset, 3)
  
  subset1 <- E.MTAB.386_eset[, id[[1]]]
  subset2 <- E.MTAB.386_eset[, id[[2]]]
  subset3 <- E.MTAB.386_eset[, id[[3]]]
  
  ## Number of observations in the set does not need to be a multiple of
  ## the fold parameter
  id2 <- cvSubsets(E.MTAB.386_eset, 5)
  subsets <- list()
  subsets[[1]] <- E.MTAB.386_eset[, id2[[1]]]
  subsets[[2]] <- E.MTAB.386_eset[, id2[[2]]]
  subsets[[3]] <- E.MTAB.386_eset[, id2[[3]]]
  subsets[[4]] <- E.MTAB.386_eset[, id2[[4]]]
  subsets[[5]] <- E.MTAB.386_eset[, id2[[5]]]
  
  ## Support matrix
  X <- exprs(subset1) 
  id3 <- cvSubsets(X, 5)
  
  ## Support SummarizedExperiment 
  nrows <- 200; ncols <- 6
  counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
  rowData <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),
                     IRanges(floor(runif(200, 1e5, 1e6)), width=100),
                     strand=sample(c("+", "-"), 200, TRUE))
  colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),
                       row.names=LETTERS[1:6])
  sset <- SummarizedExperiment(assays=SimpleList(counts=counts),
                               rowData=rowData, colData=colData)
  id4 <- cvSubsets(sset, 5)
  
})