\name{plusMinus}
\alias{plusMinus}
\title{plusMinus}
\description{function for plusMinus algorithm}
\usage{plusMinus(X, y, lambda = NULL, tuningpar = "nfeatures", standardize = FALSE, 
    directionality = "posneg", ties.method = "average", votingthresholdquantile = 0.5, 
    modeltype = "plusminus")}
\arguments{
  \item{X}{gene expression matrix}
  \item{y}{response variables}
  \item{lambda}{lambda}
  \item{tuningpar}{tuning parameter}
  \item{standardize}{standardize or not}
  \item{directionality}{directionality                       }
  \item{ties.method}{ties.method                      }
  \item{votingthresholdquantile}{votingthresholdquantile}
  \item{modeltype}{modeltype}
}

\value{returns regression coefficients }

\author{Yuqing Zhang, Christoph Bernau, Levi Waldron}




\examples{
set.seed(8)
library(curatedOvarianData)
data( E.MTAB.386_eset )
eset <- E.MTAB.386_eset[1:100, 1:30]
rm(E.MTAB.386_eset)

X <- t(exprs(eset))

time <- eset$days_to_death
cens <- sample(0:1, 30, replace=TRUE)
y <- Surv(time, cens)

beta <- plusMinus(X, y)
beta
}
