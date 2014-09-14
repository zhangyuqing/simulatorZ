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
\usage{getTrueModel(obj, y.vars, parstep)}
\arguments{
  \item{obj}{a list of ExpressionSets, matrix or SummarizedExperiment}
  \item{y.vars}{a list of response variables, Surv, matrix or data.frame object}
  \item{parstep}{number of steps in CoxBoost}
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
library(GenomicRanges)
data(GSE17260_eset)
data(E.MTAB.386_eset)
data(GSE14764_eset)
esets <- list(GSE17260=GSE17260_eset, E.MTAB.386=E.MTAB.386_eset, GSE14764=GSE14764_eset)
esets.list <- lapply(esets, function(eset){
  return(eset[1:500, 1:20])
})

## simulate on multiple ExpressionSets
set.seed(8) 

y.list <- lapply(esets.list, function(eset){
  time <- eset$days_to_death
  cens.chr <- eset$vital_status
  cens <- c()
  for(i in 1:length(cens.chr)){
    if(cens.chr[i] == "living") cens[i] <- 1
    else cens[i] <- 0
  }
  y <- Surv(time, cens)
  return(y)
})
   
res1 <- getTrueModel(esets.list, y.list, 100)
## Get true model from one set
res2 <- getTrueModel(esets.list[1], y.list[1], 100)
names(res2)
res2$lp
## note that y.list[1] cannot be replaced by y.list[[1]]

## Support matrices
X.list <- lapply(esets.list, function(eset){
  return(exprs(eset))
})
res3 <- getTrueModel(X.list, y.list, 100)  
}
