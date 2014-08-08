\name{simBootstrap}
\alias{simBootstrap}
\title{simBootstrap}
\description{the driver function to perform three-step bootstrap}
\usage{simBootstrap(esets, funSimData = simData, balance.variables = NULL, 
    n.samples, type = "simulated", funParameter = getParameter, 
    step, funSurvTime = generateSurvtime)}
\arguments{
  \item{esets}{the original ExpressionSets}
  \item{funSimData}{function to perform non-parametric bootstrap}
  \item{balance.variables}{covariate names to balance in the simulated sets}
  \item{n.samples}{number of samples to generate}
  \item{type}{whether to include parametric bootstrap}
  \item{funParameter}{function to perform CoxBoost and return cumulative hazard of original sets}
  \item{step}{step number of CoxBoost}
  \item{funSurvTime}{function to perform parametric bootstrap}
}

\value{a list of values including:
esets.list = a list of ExpressionSets, if type==simulated, then the survival time is
generated. 
indices.list = a list of indices indicating which sample the simulated sample is in the 
original set
setsID = a vector to indicate the original ID of simulated sets, if 
type=="original", setsID should be 1,2,3,...
lp.list = a list of true linear predictor of each original data sets
beta.list = a list of true coefficients used for simulating observations}

\author{Yuqing Zhang}





