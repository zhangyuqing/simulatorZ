\name{simBootstrap}
\alias{simBootstrap}
\title{simBootstrap}
\description{the driver function to perform three-step bootstrap simulation}
\usage{simBootstrap(esets, y.vars, n.samples, parstep, type = "two-steps", 
    balance.variables = NULL, funSimData = simData, funTrueModel = getTrueModel, 
    funSurvTime = simTime)}
\arguments{
  \item{esets}{a list of ExpressionSets }
  \item{y.vars}{a list of reponse variables}
  \item{n.samples}{number of samples to resample in each set}
  \item{parstep}{step number to fit CoxBoost}
  \item{type}{whether to include resampling set labels }
  \item{balance.variables}{covariate names to balance in the simulated sets}
  \item{funSimData}{function to perform non-parametric bootstrap}
  \item{funTrueModel}{function to construct true models in original sets}
  \item{funSurvTime}{function to perform parametric bootstrap}
}

\value{a list of values including:
esets.list = a list of simulated ExpressionSets
indices.list = a list of indices indicating which sample the simulated sample is in the 
original set
setsID = a vector to indicate the original ID of simulated sets, if 
type=="original", setsID should be 1,2,3,...
lp.list = a list of true linear predictor of each original data sets
beta.list = a list of true coefficients used for simulating observations
survH.list = list of cumulative survival hazard
censH.list = list of cumulative censoring hazard
grid.list = list of timeline grid corresponding to survH and censH respectivley}

\author{Levi Waldron, Christoph Bernau, Yuqing Zhang}




\examples{
library(curatedOvarianData)
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

simmodels <- simBootstrap(esets=esets, y.vars=y.vars, 10, 100)

# skip resampling labels
simmodels <- simBootstrap(esets=esets, y.vars=y.vars, 10, 100,
                          type="one-step")

# balance covariates
simmodels <- simBootstrap(esets=esets, y.vars=y.vars, 10, 100,
                          balance.variables="tumorstage")
}
