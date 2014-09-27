\name{cvSubsets}
\alias{cvSubsets}
\title{cvSubsets}
\description{To generate a list of subsets(indices of observations) from one set}
\usage{cvSubsets(obj, fold)}
\arguments{
  \item{obj}{a ExpressionSet, matrix or SummarizedExperiment object. If it is a matrix,
columns represent samples}
  \item{fold}{the number of folds in cross validation. 
Number of observations in the set does not need to be a multiple of fold}
}

\value{returns the list of indices of subsets}

\author{Yuqing Zhang}




\examples{
library(curatedOvarianData)
data(E.MTAB.386_eset)

set.seed(8)
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
}
