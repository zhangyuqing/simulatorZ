\name{cvSubsets}
\alias{cvSubsets}
\title{cvSubsets}
\description{To generate a list of subsets(indices of observations) from one set}
\usage{cvSubsets(obj, fold)}
\arguments{
  \item{obj}{a ExpressionSet, matrix or RangedSummarizedExperiment object. If it is a
matrix,columns represent samples}
  \item{fold}{the number of folds in cross validation. 
Number of observations in the set does not need to be a multiple of fold}
}

\value{returns the list of indices of subsets}

\author{Yuqing Zhang, Christoph Bernau, Levi Waldron}




\examples{
library(curatedOvarianData)
data(E.MTAB.386_eset)

id <- cvSubsets(E.MTAB.386_eset, 3)
subsets <- lapply(1:3, function(i){E.MTAB.386_eset[1:10, id[[i]]]})
sapply(subsets, dim)
rm(subsets)

## Number of observations in the set does not need to be a multiple of
## the fold parameter
id2 <- cvSubsets(E.MTAB.386_eset, 5)
subsets <- lapply(1:5, function(j){E.MTAB.386_eset[1:10, id2[[j]]]})
sapply(subsets, dim)
rm(subsets)
}
