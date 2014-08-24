\name{cvSubsets}
\alias{cvSubsets}
\title{cvSubsets}
\description{To generate a list of subsets(indices of observations) from one set}
\usage{cvSubsets(eset, fold)}
\arguments{
  \item{eset}{a ExpressionSet, matrix or SummarizedExperiment object. If it is a matrix,
columns represent samples}
  \item{fold}{the number of folds in cross validation. 
Number of observations in the set does not need to be a multiple of fold}
}

\value{returns the list of indices of subsets}

\author{Levi Waldron, Christoph Bernau, Yuqing Zhang}




\examples{
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

}
