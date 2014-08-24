\name{simTime}
\alias{simTime}
\title{simTime}
\description{simTime is a function to perform the parametric-bootstrap step, where we use the true coefficients
and cumulative hazard to simulate survival and censoring.}
\usage{simTime(simmodels, result)}
\arguments{
  \item{simmodels}{a list in the form of the return value of simData()
which consists of three lists:
esets: a list of ExpressionSets, matrics or SummarizedExperiments
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

\author{Levi Waldron, Christoph Bernau, Yuqing Zhang}




\examples{
library(curatedOvarianData)
library(GenomicRanges)
data( E.MTAB.386_eset )
eset1 <- E.MTAB.386_eset[1:10, 1:5]
eset2 <- E.MTAB.386_eset[1:10, 6:10]
eset3 <- E.MTAB.386_eset[1:10, 11:15]  
esets <- list(eset1, eset2, eset3) 

time1 <- eset1$days_to_death
cens1 <- c(0, 0, 0, 1, 1)
y1 <- Surv(time1, cens1)
time2 <- eset2$days_to_death
cens2 <- c(1, 1, 0, 0, 0)
y2 <- Surv(time2, cens2)
time3 <- eset3$days_to_death
cens3 <- c(1, 0, 0, 0, 1)
y3 <- Surv(time3, cens3)
y.vars <- list(y1, y2, y3)

# To perform both parametric and non-parametric bootstrap, you can call simBootstrap()
# or, you can divide the steps into:
res <- getTrueModel(esets, y.vars, 100)
simmodels <- simData(esets=esets, y.vars=y.vars, n.samples=10)

# Then, use this function
simmodels <- simTime(simmodels=simmodels, result=res) 

# it also supports performing only the parametrc bootstrap step on a list of expressionsets
# but you need to construct the parameter by scratch
res <- getTrueModel(esets, y.vars, 100)
setsID <- 1:length(esets)
indices <- list()
for(i in setsID){
  indices[[i]] <- 1:length(sampleNames(esets[[i]])) 
}
simmodels <- list(esets=esets, y.vars=y.vars, indices=indices, setsID=setsID)

new.simmodels <- simTime(simmodels=simmodels, result=res)  
}
