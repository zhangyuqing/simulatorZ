geneFilter <- structure(function
### the function to filter genes by intergrative correlation
(esets,
 ### a list of ExpressionSet, matrix or SummarizedExperiment objects. If its elements are matrices,
 ### columns represent samples, rows represent genes
 cor.cutoff=0.5 
 ### the cutoff threshold for filtering genes
){
  geneid <- c()
  if(class(esets[[1]])=="ExpressionSet"){
    X.list <- lapply(esets, function(eset){
      return(t(exprs(eset)))
    })
  }    
  else if(class(esets[[1]])=="matrix"){
    X.list <- lapply(esets, function(eset){
      return(t(eset))
    })
  }
  else if(class(esets[[1]])=="SummarizedExperiment"){
    X.list <- lapply(esets, function(eset){
      return(t(assay(eset)))
    })
  }

  index <- 1
  qua.id <- list()
  for(i in 1:length(esets)){
    for(j in i:length(esets)){
      cor.mi <- cor(X.list[[i]])
      cor.mj <- cor(X.list[[j]])
      int.score <- diag(cor(cor.mi, cor.mj))
      qua.id[[index]] <- as.numeric(which(int.score > cor.cutoff))
      index <- index + 1
    }
  }
  
  geneid <- 1:ncol(X.list[[1]])
  for(k in 1:length(qua.id)){
    geneid <- intersect(geneid, qua.id[[k]])
  }
    
  new.esets <- lapply(esets, function(eset){
    return(eset[geneid,])
  })
  
  return(new.esets)
  ### returns a list of ExpressionSets with genes filtered 
},ex=function(){
  library(curatedOvarianData)
  library(GenomicRanges)
  data( E.MTAB.386_eset )
  eset1 <- E.MTAB.386_eset[1:500, 1:5]
  eset2 <- E.MTAB.386_eset[1:500, 6:10]
  eset3 <- E.MTAB.386_eset[1:500, 11:15]  
  esets <- list(eset1, eset2, eset3) 
  
  result.set <- geneFilter(esets, 0)
  result.set
  ### as we cannot calculate correlation with one set, this function just 
  ### delivers the same set if esets has length 1
  result.oneset <- geneFilter(list(eset1))
  result.oneset
  
  ## Support matrices
  X.list <- lapply(esets, function(eset){
    return(exprs(eset)) ## Columns represent samples!
  })
  result.set <- geneFilter(X.list, 0)
  dim(result.set[[1]])
  
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
  s.list <- list(sset, sset)
  result.set <- geneFilter(s.list, 0.9) 
  ## the same set should resemble each other, no genes filtered
  dim(assay(result.set[[1]]))
})
