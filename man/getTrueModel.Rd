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
\usage{getTrueModel(esets, y.vars, parstep, balance.variables = NULL)}
\arguments{
  \item{esets}{a list of ExpressionSets, matrix or SummarizedExperiment}
  \item{y.vars}{a list of response variables}
  \item{parstep}{CoxBoost parameter}
  \item{balance.variables}{variable names to be balanced.}
}

\value{returns a list of values:
beta: True coefficients obtained by fitting CoxBoost to the original ExpressionSets
grid: timeline grid corresponding to hazard estimations censH and survH
survH: cumulative hazard for survival times distribution
censH: cumulative hazard for censoring times distribution
lp: true linear predictors }

\author{Yuqing Zhang, Christoph Bernau, Levi Waldron}




\examples{
library(curatedOvarianData)
data(GSE14764_eset)
data(E.MTAB.386_eset)
esets.list <- list(GSE14764=GSE14764_eset[1:500, 1:20], 
                   E.MTAB.386=E.MTAB.386_eset[1:500, 1:20])
rm(E.MTAB.386_eset, GSE14764_eset)

## simulate on multiple ExpressionSets
set.seed(8) 

y.list <- lapply(esets.list, function(eset){
  time <- eset$days_to_death
  cens.chr <- eset$vital_status
  cens <- rep(0, length(cens.chr))
  cens[cens.chr=="living"] <- 1
  return(Surv(time, cens))
})
   
res1 <- getTrueModel(esets.list, y.list, 100)
## Get true model from one set
res2 <- getTrueModel(esets.list[1], y.list[1], 100)
names(res2)
res2$lp
## note that y.list[1] cannot be replaced by y.list[[1]]
}
