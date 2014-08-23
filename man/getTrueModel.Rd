\name{getTrueModel}
\alias{getTrueModel}
\title{getTrueModel}
\description{The parametric bootstrap simulation depends on the true model of original sets.
This function is to generate useful values from the true models for further analysis.  
We fit CoxBoost to the original sets and use the coefficients to simulate
the survival and censoring time. grid, survH, censH, which are useful for this purpose.
grid=grid corresponding to hazard estimations censH and survH
survH=cumulative hazard for survival times distribution
censH=cumulative hazard for censoring times distribution}
\usage{getTrueModel(esets, y.vars, parstep)}
\arguments{
  \item{esets}{a list of ExpressionSets}
  \item{y.vars}{a list of response variables}
  \item{parstep}{CoxBoost parameter}
}

\value{returns a list of values:
beta: True coefficients obtained by fitting CoxBoost to the original ExpressionSets
grid: timeline grid corresponding to hazard estimations censH and survH
survH: cumulative hazard for survival times distribution
censH: cumulative hazard for censoring times distribution
lp: true linear predictors }

\author{Levi Waldron, Christoph Bernau, Yuqing Zhang}




\examples{
library(curatedOvarianData)
data( E.MTAB.386_eset )
eset1 <- E.MTAB.386_eset[, 1:5]
eset2 <- E.MTAB.386_eset[, 6:10]
eset3 <- E.MTAB.386_eset[, 11:15]

## simulate on multiple ExpressionSets
esets.list <- list(eset1, eset2, eset3) 

time1 <- eset1$days_to_death
cens1 <- c(0, 0, 0, 1, 1)
y1 <- Surv(time1, cens1)
time2 <- eset2$days_to_death
cens2 <- c(1, 1, 0, 0, 0)
y2 <- Surv(time2, cens2)
time3 <- eset3$days_to_death
cens3 <- c(1, 0, 0, 0, 1)
y3 <- Surv(time3, cens3)
y.list<- list(y1, y2, y3) 

 
res1 <- getTrueModel(esets.list, y.list, 100)
res2 <- getTrueModel(list(eset1), y.list[1], 100)
names(res2)
res2$lp
## note that y.list[1] cannot be replaced by y.list[[1]]

}
