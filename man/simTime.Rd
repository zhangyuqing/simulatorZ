\name{simTime}
\alias{simTime}
\title{simTime}
\description{simTime is a function to perform the parametric-bootstrap step, where we use the true coefficients
and cumulative hazard to simulate survival and censoring.}
\usage{simTime(simmodels, result)}
\arguments{
  \item{simmodels}{a list in the form of the return value of simData()
which consists of three lists:
obj: a list of ExpressionSets, matrics or SummarizedExperiments
setsID: a list of set labels indicating which original set the simulated one is from
indices: a list of patient labels to tell which patient in the original set is drawn}
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
library(GenomicRanges)
source(system.file("extdata", "patientselection.config",
                   package="curatedOvarianData"))
source(system.file("extdata", "createEsetList.R", package="curatedOvarianData"))
esets.list <- lapply(esets, function(eset){
  return(eset[1:500, 1:20])
})
esets.list <- esets.list[1:5]
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

# To perform both parametric and non-parametric bootstrap, you can call simBootstrap()
# or, you can divide the steps into:
res <- getTrueModel(esets.list, y.list, 100)
simmodels <- simData(obj=esets.list, y.vars=y.list, n.samples=10)

# Then, use this function
simmodels <- simTime(simmodels=simmodels, result=res) 

# it also supports performing only the parametrc bootstrap step on a list of expressionsets
# but you need to construct the parameter by scratch
res <- getTrueModel(esets.list, y.list, 100)
setsID <- 1:length(esets.list)
indices <- list()
for(i in setsID){
  indices[[i]] <- 1:length(sampleNames(esets.list[[i]])) 
}
simmodels <- list(obj=esets.list, y.vars=y.list, indices=indices, setsID=setsID)

new.simmodels <- simTime(simmodels=simmodels, result=res)  
}
