\name{funCV}
\alias{funCV}
\title{funCV}
\description{Cross validation function}
\usage{funCV(eset, fold, y.var, trainFun = masomenos, funCvSubset = cvSubsets)}
\arguments{
  \item{eset}{one ExpressionSet to do cross validation with }
  \item{fold}{the number of folds in cross validation}
  \item{y.var}{response variable, matrix, data.frame(with 2 columns) or Surv object}
  \item{trainFun}{training function, which takes gene expression matrix X and response variable y as input, the coefficients as output }
  \item{funCvSubset}{function to divide one Expression Set into subsets for cross validation}
}

\value{returns the c statistics of cross validation(CV)}

\author{Levi Waldron, Christoph Bernau, Yuqing Zhang}




\examples{
library(curatedOvarianData)
data( E.MTAB.386_eset )
eset <- E.MTAB.386_eset[1:100, 1:30]
time <- sample(2000:5000, 30)
cens <- sample(0:1, 30, replace=TRUE)
y <- Surv(time, cens)  
y1 <- cbind(time, cens)

funCV(eset, 3, y)
funCV(eset, 3, y1, trainFun=plusMinus)
## any training function will do as long as it takes the gene expression matrix X
## and response variable y(matrix, data.frame or Surv object) as parameters, and
## return the coefficients as its value
}
