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
  \item{lambda}{
}
  \item{tuningpar}{
}
  \item{standardize}{
}
  \item{directionality}{
}
  \item{ties.method}{
}
  \item{votingthresholdquantile}{
}
  \item{modeltype}{
}
}

\value{returns regression coefficients }

\author{Levi Waldron, Christoph Bernau, Yuqing Zhang}




\examples{
library(curatedOvarianData)
data( E.MTAB.386_eset )
eset <- E.MTAB.386_eset[1:100, 1:30]

X <- t(exprs(eset))

time <- eset$days_to_death
cens <- sample(0:1, 30, replace=TRUE)
y <- Surv(time, cens)

beta <- plusMinus(X, y)
beta
}
