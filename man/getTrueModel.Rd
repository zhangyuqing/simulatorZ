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

y.list <- list()
for(i in 1:length(esets.list)){
  time <- esets.list[[i]]$days_to_death
  cens <- sample(0:1, 5, replace=TRUE)
  y.list[[i]] <- Surv(time, cens)
}  

res1 <- getTrueModel(esets.list, y.list, 100)
res2 <- getTrueModel(list(eset1), y.list[1], 100)
## note that y.list[1] cannot be replaced by y.list[[1]]

}
