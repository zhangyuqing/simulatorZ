\name{simTime}
\alias{simTime}
\title{simTime}
\description{simTime is a function to perform the parametric-bootstrap step, where we use the true coefficients
and cumulative hazard to simulate survival and censoring.}
\usage{simTime(simmodels, original.yvars, result)}
\arguments{
  \item{simmodels}{a list in the form of the return value of simData()
which consists of three lists:
obj: a list of ExpressionSets, matrices or RangedSummarizedExperiments
setsID: a list of set labels indicating which original set the simulated one is from
indices: a list of patient labels to tell which patient in the original set is drawn}
  \item{original.yvars}{response variable in the order of original sets(without sampling)}
  \item{result}{a list in the form of return of getTrueModel()
which consists of five lists: 
Beta: a list of coefficients obtained by 
grid: timeline grid corresponding to hazard estimations censH and survH
survH: cumulative hazard for survival times distribution
censH: cumulative hazard for censoring times distribution
lp: true linear predictors}
}

\value{survival time is saved in phenodata, here the function still returns the ExpressionSets}

\author{Yuqing Zhang, Christoph Bernau, Levi Waldron}




\examples{
library(curatedOvarianData)
data(E.MTAB.386_eset)
data(GSE14764_eset)
esets.list <- list(E.MTAB.386=E.MTAB.386_eset[1:100, 1:20], GSE14764=GSE14764_eset[1:100, 1:20])
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

# To perform both parametric and non-parametric bootstrap, you can call simBootstrap()
# or, you can divide the steps into:
res <- getTrueModel(esets.list, y.list, 100)
simmodels <- simData(obj=esets.list, y.vars=y.list, n.samples=10)

# Then, use this function
simmodels <- simTime(simmodels=simmodels, original.yvars=y.list, result=res) 

# it also supports performing only the parametrc bootstrap step on a list of expressionsets
# but you need to construct the parameter by scratch
res <- getTrueModel(esets.list, y.list, 100)
setsID <- seq_along(esets.list)
indices <- list()
for(i in setsID){
  indices[[i]] <- seq_along(sampleNames(esets.list[[i]])) 
}
simmodels <- list(obj=esets.list, y.vars=y.list, indices=indices, setsID=setsID)

new.simmodels <- simTime(simmodels=simmodels, original.yvars=y.list, result=res)  
}
