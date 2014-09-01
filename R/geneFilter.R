geneFilter <- structure(function
### the function to filter genes by intergrative correlation
(esets,
 ### a list of ExpressionSet, matrix or SummarizedExperiment objects. If its elements are matrices,
 ### columns represent samples, rows represent genes
 cor.cutoff=0.5 
 ### the cutoff threshold for filtering genes
){
  index <- 1
  qua.id <- list()
  geneid <- c()
  
  for(i in 1:length(esets)){
    for(j in i:length(esets)){
      if(i!=j){
        print(paste(i, j, sep="--"))
        if(class(esets[[1]])=="ExpressionSet"){
          Xi <- t(exprs(esets[[i]]))
          Xj <- t(exprs(esets[[j]]))
        }    
        else if(class(esets[[1]])=="matrix"){
          Xi <- t(esets[[i]])
          Xj <- t(esets[[j]])
        }
        else if(class(esets[[1]])=="SummarizedExperiment"){
          Xi <- t(assay(esets[[i]]))
          Xj <- t(assay(esets[[j]]))
        }   
        m1 <- cor(Xi)
        m2 <- cor(Xj)
        int.score <- c()
        for(k in 1:ncol(m1)){
          int.score[k] <- cor(m1[, k], m2[, k])
        }
        qua.id[[index]] <- as.numeric(which(int.score > cor.cutoff))
        index <- index + 1
        geneid <- 1:ncol(Xi)
        rm(Xi, Xj, m1, m2)
      }      
    }
  }
    
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
