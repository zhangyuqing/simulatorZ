getMatrix <- structure(function
### a generic function with methods for matrix, ExpressionSet, and RangedSummarizedExperiment objects.
(obj
### Matrix, ExpressionSet or RangedSummarizedExperiment
){
    if (is.matrix(obj))
        return(obj)
    if (is(obj, "ExpressionSet"))
        return(exprs(obj))
    if (is(obj, "RangedSummarizedExperiment")
     || is(obj, "SummarizedExperiment"))  # for backward compatibility
        return(assay(obj))
	### return the expression matrix
    stop("Wrong class of obj!")
},ex=function(){
  library(curatedOvarianData)
  data( E.MTAB.386_eset )
  eset <- E.MTAB.386_eset[1:100, 1:30]
  
  X <- getMatrix(eset)
})

