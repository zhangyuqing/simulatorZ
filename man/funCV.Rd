\name{funCV}
\alias{funCV}
\title{funCV}
\description{Cross validation function}
\usage{funCV(obj, fold, y.var, trainFun = masomenos, funCvSubset = cvSubsets)}
\arguments{
  \item{obj}{a ExpressionSet, matrix or SummarizedExperiment object. If it is a matrix,
columns represent samples }
  \item{fold}{the number of folds in cross validation}
  \item{y.var}{response variable, matrix, data.frame(with 2 columns) or Surv object}
  \item{trainFun}{training function, which takes gene expression matrix X and response variable y as input, the coefficients as output }
  \item{funCvSubset}{function to divide one Expression Set into subsets for cross validation}
}

\value{returns the c statistics of cross validation(CV)}

\author{Yuqing Zhang}




\examples{
library(curatedOvarianData)
library(GenomicRanges)
set.seed(8)
data( E.MTAB.386_eset )
eset <- E.MTAB.386_eset[1:100, 1:30]
time <- eset$days_to_death
cens.chr <- eset$vital_status
cens <- c()
for(i in 1:length(cens.chr)){
  if(cens.chr=="living") cens[i] <- 1
  else cens[i] <- 0
}
y <- Surv(time, cens)  
y1 <- cbind(time, cens)

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

funCV(eset, 3, y)
funCV(eset, 3, y1, trainFun=plusMinus)
funCV(exprs(eset), 3, y)

funCV(sset, 3, y.vars)
## any training function will do as long as it takes the gene expression matrix X
## and response variable y(matrix, data.frame or Surv object) as parameters, and
## return the coefficients as its value
}
