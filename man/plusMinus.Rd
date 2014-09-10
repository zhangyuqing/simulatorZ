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
\references{Zhao, S., Huttenhower, G. P. C., and Waldron, L. (2013). Mas-o-menos:
a simple sign averaging method for discrimination in genomic data analysis.
http://biostats.bepress.com/harvardbiostat/paper158/. Accessed: 2013-10-24.}
\author{Yuqing Zhang, Christoph Bernau, Levi Waldron}




\examples{
set.seed(8)
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
