## Alternatively, this could be implemented as a generic function with methods
## for matrix, ExpressionSet, and RangedSummarizedExperiment objects.
getMatrix <- function(obj)
{
    if (is.matrix(obj))
        return(obj)
    if (is(obj, "ExpressionSet"))
        return(exprs(obj))
    if (is(obj, "RangedSummarizedExperiment")
     || is(obj, "SummarizedExperiment"))  # for backward compatibility
        return(assay(obj))
    stop("Wrong class of obj!")
}

