\name{geneFilter}
\alias{geneFilter}
\title{geneFilter}
\description{the function to filter genes by intergrative correlation}
\usage{geneFilter(esets, cor.cutoff = 0.6)}
\arguments{
  \item{esets}{a list of Expression Sets with same number of observations and genes }
  \item{cor.cutoff}{the cutoff threshold for filtering genes}
}

\value{returns a list of ExpressionSets with genes filtered }

\author{Levi Waldron, Christoph Bernau, Yuqing Zhang}




\examples{
library(curatedOvarianData)
data( E.MTAB.386_eset )
eset1 <- E.MTAB.386_eset[, 1:5]
eset2 <- E.MTAB.386_eset[, 6:10]
eset3 <- E.MTAB.386_eset[, 11:15]  
esets <- list(eset1, eset2, eset3) 

result.set <- geneFilter(esets)
result.set
### as we cannot calculate correlation with one set, this function just 
### delivers the same set if esets has length 1
result.oneset <- geneFilter(list(eset1))
result.oneset
}
