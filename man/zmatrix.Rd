\name{zmatrix}
\alias{zmatrix}
\title{zmatrix}
\description{generate a matrix of c statistics}
\usage{zmatrix(esets, y.vars, fold, trainingFun = masomenos, cvFun = funCV, 
    cvSubsetFun = cvSubsets)}
\arguments{
  \item{esets}{a list of ExpressionSet, matrix or SummarizedExperiment objects, if its 
elements are matrices, columns represent samples }
  \item{y.vars}{a list of response variables, all the response variables shold be
matrix, data.frame(with 2 columns) or Surv object}
  \item{fold}{cvFun parameter, in this case passes to funCV()}
  \item{trainingFun}{training function}
  \item{cvFun}{function to perform cross study within one set}
  \item{cvSubsetFun}{function to divide the expression sets into subsets for cross validation }
}

\value{outputs one matrix of validation statistics}

\author{Levi Waldron, Christoph Bernau, Yuqing Zhang}




\examples{
library(curatedOvarianData)
library(GenomicRanges)
data( E.MTAB.386_eset )
eset1 <- E.MTAB.386_eset[1:100, 1:30]
eset2 <- E.MTAB.386_eset[1:100, 31:60]
eset3 <- E.MTAB.386_eset[1:100, 61:90]  
esets <- list(eset1, eset2, eset3) 

time1 <- eset1$days_to_death
#cens1 <- c(0, 0, 0, 1, 1)
cens1 <- sample(0:1, 30, replace=TRUE)
y1 <- Surv(time1, cens1)
time2 <- eset2$days_to_death
#cens2 <- c(1, 1, 0, 0, 0)
cens2 <- sample(0:1, 30, replace=TRUE)
y2 <- Surv(time2, cens2)
time3 <- eset3$days_to_death
#cens3 <- c(1, 0, 0, 0, 1)
cens3 <- sample(0:1, 30, replace=TRUE)
y3 <- Surv(time3, cens3)
y.vars <- list(y1, y2, y3)

# generate on original ExpressionSets
z <- zmatrix(esets, y.vars, 3)

# generate on simulated ExpressionSets
simmodels <- simBootstrap(esets, y.vars, 10, 100)
z <- zmatrix(simmodels$esets.list, simmodels$y.vars.list, 3)

# support matrix
X.list <- lapply(esets, function(eset){
  newx <- exprs(eset) ### columns represent samples !!
  return(newx)
}) 
z <- zmatrix(X.list, y.vars, 3)

# support SummarizedExperiment
nrows <- 200; ncols <- 6
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
rowData <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),
                   IRanges(floor(runif(200, 1e5, 1e6)), width=100),
                   strand=sample(c("+", "-"), 200, TRUE))
colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),
                     row.names=LETTERS[1:6])
sset <- SummarizedExperiment(assays=SimpleList(counts=counts),
                             rowData=rowData, colData=colData)

time <- sample(4500:4700, 6, replace=TRUE)
cens <- sample(0:1, 6, replace=TRUE)
y.vars <- Surv(time, cens)

z <- zmatrix(list(sset[,1:3], sset[,4:6]), list(y.vars[1:3,],y.vars[4:6,]), 3)
}
