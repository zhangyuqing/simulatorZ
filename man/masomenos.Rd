\name{masomenos}
\alias{masomenos}
\title{masomenos}
\description{function for Mas-o-menos algorithm}
\usage{masomenos(X, y, option = "fast", ...)}
\arguments{
  \item{X}{matrix with rows corresponding to subjects and columns to features resp}
  \item{y}{response variable, a data.frame, matrix, or Surv object: c(time, event)}
  \item{option}{whether to use C or R code to fit the marginal Cox models}
  \item{\dots}{
}
}

\value{return the coefficients}
\references{Zhao, S., Huttenhower, G. P. C., and Waldron, L. (2013). Mas-o-menos:
a simple sign averaging method for discrimination in genomic data analysis.
http://biostats.bepress.com/harvardbiostat/paper158/. Accessed: 2013-10-24.}
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

beta <- masomenos(X=X, y=y)
beta
}
