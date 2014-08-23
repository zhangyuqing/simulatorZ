geneFilter <- structure(function
### the function to filter genes by intergrative correlation
(esets,
 ### a list of Expression Sets with same number of observations and genes 
 cor.cutoff=0.5 
 ### the cutoff threshold for filtering genes
){
  geneid <- c()
  X.list <- lapply(esets, function(eset){
    return(t(exprs(eset)))
  }) 

  index <- 1
  for(k in 1:ncol(X.list[[1]])){
    marker <- TRUE
    for(i in 1:length(esets)){
      for(j in i:length(esets)){
        cor.mi <- cor(X.list[[i]])
        cor.mj <- cor(X.list[[j]])
        if(cor(cor.mi[k, ], cor.mj[k, ]) < cor.cutoff) marker <- FALSE
        rm(cor.mi, cor.mj)
      }
    }
    if(marker) {
      geneid[index] <- k
      index <- index + 1
    }
  }
  
  
  new.esets <- lapply(esets, function(eset){
    return(eset[geneid,])
  })
  
  return(new.esets)
  ### returns a list of ExpressionSets with genes filtered 
},ex=function(){
  library(curatedOvarianData)
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
})
