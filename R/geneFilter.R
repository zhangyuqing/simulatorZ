geneFilter <- structure(function
### the function to filter genes by Intergrative Correlation
##references<< Garrett-Mayer, E., Parmigiani, G., Zhong, X., Cope, L., 
## Gabrielson, E., Cross-study validation and combined analysis of gene 
## expression microarray data. Biostatistics. 2008 Apr;9(2):333-354.
(obj,
 ### a list of ExpressionSet, matrix or SummarizedExperiment objects. If its elements are matrices,
 ### columns represent samples, rows represent genes
 cor.cutoff=0.5 
 ### the cutoff threshold for filtering genes. Only when the integrative correlation
 ### between every pair of sets is larger than the cutoff value, will the gene 
 ### be selected.
){
  index <- 1
  qua.id <- list()
  geneid <- c()
  
  for(i in 1:length(obj)){
    for(j in i:length(obj)){
      if(i!=j){
        print(paste(i, j, sep="--"))
        if(class(obj[[1]])=="ExpressionSet"){
          Xi <- t(exprs(obj[[i]]))
          Xj <- t(exprs(obj[[j]]))
        }    
        else if(class(obj[[1]])=="matrix"){
          Xi <- t(obj[[i]])
          Xj <- t(obj[[j]])
        }
        else if(class(obj[[1]])=="SummarizedExperiment"){
          Xi <- t(assay(obj[[i]]))
          Xj <- t(assay(obj[[j]]))
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
  
  if(length(obj)==1){
    if(class(obj[[1]])=="ExpressionSet"){
      X <- t(exprs(obj[[1]]))
      X <- t(exprs(obj[[1]]))
    }    
    else if(class(obj[[1]])=="matrix"){
      X <- t(obj[[1]])
      X <- t(obj[[1]])
    }
    else if(class(obj[[1]])=="SummarizedExperiment"){
      X <- t(assay(obj[[1]]))
      X <- t(assay(obj[[1]]))
    }
    geneid <- qua.id[[1]] <- 1:ncol(X)
    rm(X)
  }
  
  for(k in 1:length(qua.id)){
    geneid <- intersect(geneid, qua.id[[k]])
  }
    
  new.obj <- lapply(obj, function(obj.ele){
    return(obj.ele[geneid,])
  })
  
  return(new.obj)
  ### returns a list of ExpressionSets matrix or SummarizedExperiment objects
  ### with genes filtered 
},ex=function(){
  set.seed(8)
  library(curatedOvarianData)
  library(GenomicRanges)
  data(GSE17260_eset)
  data(E.MTAB.386_eset)
  data(GSE14764_eset)
  esets <- list(GSE17260=GSE17260_eset, E.MTAB.386=E.MTAB.386_eset, GSE14764=GSE14764_eset)
  esets.list <- lapply(esets, function(eset){
    return(eset[1:1500, 1:10])
  })
  
  result.set <- geneFilter(esets.list, 0)
  result.set
  ### as we cannot calculate correlation with one set, this function just 
  ### delivers the same set if esets has length 1
  result.oneset <- geneFilter(esets.list[1])
  result.oneset
  
  ## Support matrices
  X.list <- lapply(esets.list, function(eset){
    return(exprs(eset)) ## Columns represent samples!
  })
  result.set <- geneFilter(X.list, 0)
  dim(result.set[[1]])
  
  ## Support SummarizedExperiment
  nrows <- 200; ncols <- 6
  counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
  rowRanges <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),
                       IRanges(floor(runif(200, 1e5, 1e6)), width=100),
                       strand=sample(c("+", "-"), 200, TRUE))
  colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),
                       row.names=LETTERS[1:6])
  sset <- SummarizedExperiment(assays=SimpleList(counts=counts),
                               rowRanges=rowRanges, colData=colData)
  s.list <- list(sset, sset)
  result.set <- geneFilter(s.list, 0.9) 
  ## the same set should resemble each other, no genes filtered
  dim(assay(result.set[[1]]))
})
