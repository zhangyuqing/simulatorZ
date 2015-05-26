### a generic function with methods for matrix, ExpressionSet, and RangedSummarizedExperiment objects.
### Matrix, ExpressionSet or RangedSummarizedExperiment

getMatrix <- function(obj){
    if (is.matrix(obj))
        return(obj)
    if (is(obj, "ExpressionSet"))
        return(exprs(obj))
    if (is(obj, "RangedSummarizedExperiment")
     || is(obj, "SummarizedExperiment"))  # for backward compatibility
        return(assay(obj))
	### return the expression matrix
    stop("Wrong class of obj!")
}

