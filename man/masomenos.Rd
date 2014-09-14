\name{masomenos}
\alias{masomenos}
\title{masomenos}
\description{function for Mas-o-menos algorithm}
\usage{masomenos(X, y)}
\arguments{
  \item{X}{gene expression matrix, with rows corresponding to patients and columns corresponding to genes }
  \item{y}{response variable, a data.frame, matrix or Surv object, with the first column as time to event and the second column 
censoring status }
}

\value{return betas}
\references{Zhao, S., Huttenhower, G. P. C., and Waldron, L. (2013). Mas-o-menos:
a simple sign averaging method for discrimination in genomic data analysis.
http://biostats.bepress.com/harvardbiostat/paper158/. Accessed: 2013-10-24.}
\author{Yuqing Zhang}




\examples{
set.seed(8)
library(curatedOvarianData)
data( E.MTAB.386_eset )
eset <- E.MTAB.386_eset[1:100, 1:30]

X <- t(exprs(eset))

time <- eset$days_to_death
cens <- sample(0:1, 30, replace=TRUE)
y <- Surv(time, cens)

beta <- masomenos(X, y)
beta
}
