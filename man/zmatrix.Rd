\name{zmatrix}
\alias{zmatrix}
\title{zmatrix}
\description{generate a matrix of c statistics}
\usage{zmatrix(obj, y.vars, fold, trainingFun = masomenos, cvFun = funCV, 
    cvSubsetFun = cvSubsets)}
\arguments{
  \item{obj}{a list of ExpressionSet, matrix or RangedSummarizedExperiment objects, if its 
elements are matrices, columns represent samples }
  \item{y.vars}{a list of response variables, all the response variables shold be
matrix, data.frame(with 2 columns) or Surv object}
  \item{fold}{cvFun parameter, in this case passes to funCV()}
  \item{trainingFun}{training function}
  \item{cvFun}{function to perform cross study within one set}
  \item{cvSubsetFun}{function to divide the expression sets into subsets for cross validation }
}

\value{outputs one matrix of validation statistics}

\author{Yuqing Zhang, Christoph Bernau, Levi Waldron}




\examples{
library(curatedOvarianData)
library(GenomicRanges)
data(GSE17260_eset)
data(E.MTAB.386_eset)
data(GSE14764_eset)
esets <- list(GSE17260=GSE17260_eset, E.MTAB.386=E.MTAB.386_eset, GSE14764=GSE14764_eset)
esets.list <- lapply(esets, function(eset){
  return(eset[1:500, 1:20])
})

## simulate on multiple ExpressionSets
set.seed(8) 

y.list <- lapply(esets.list, function(eset){
  time <- eset$days_to_death
  cens.chr <- eset$vital_status
  cens <- c()
  for(i in seq_along(cens.chr)){
    if(cens.chr[i] == "living") cens[i] <- 1
    else cens[i] <- 0
  }
  y <- Surv(time, cens)
  return(y)
})

# generate on original ExpressionSets
z <- zmatrix(esets.list, y.list, 3)

# generate on simulated ExpressionSets
simmodels <- simBootstrap(esets.list, y.list, 10, 100)
z <- zmatrix(simmodels$obj.list, simmodels$y.vars.list, 3)

# support matrix
X.list <- lapply(esets.list, function(eset){
  newx <- exprs(eset) ### columns represent samples !!
  return(newx)
}) 
z <- zmatrix(X.list, y.list, 3)

# support RangedSummarizedExperiment
nrows <- 200; ncols <- 6
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
rowRanges <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),
                     IRanges(floor(runif(200, 1e5, 1e6)), width=100),
                     strand=sample(c("+", "-"), 200, TRUE))
colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),
                     row.names=LETTERS[1:6])
sset <- SummarizedExperiment(assays=SimpleList(counts=counts),
                             rowRanges=rowRanges, colData=colData)

time <- sample(4500:4700, 6, replace=TRUE)
cens <- sample(0:1, 6, replace=TRUE)
y.vars <- Surv(time, cens)

z <- zmatrix(list(sset[,1:3], sset[,4:6]), list(y.vars[1:3,],y.vars[4:6,]), 3)
}
