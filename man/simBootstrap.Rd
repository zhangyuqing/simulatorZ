\name{simBootstrap}
\alias{simBootstrap}
\title{simBootstrap}
\description{the driver function to perform three-step bootstrap}
\usage{simBootstrap(esets, funSimData, balance.variables, funParameter, 
    par_step, funSurvTime)}
\arguments{
  \item{esets}{the original ExpressionSets}
  \item{funSimData}{function to perform non-parametric bootstrap}
  \item{balance.variables}{covariate names to balance in the simulated sets}
  \item{funParameter}{function to perform CoxBoost and return cumulative hazard of original sets}
  \item{par_step}{step number of CoxBoost}
  \item{funSurvTime}{function to perform parametric bootstrap}
}

\value{a list of simulated ExpressionSets}

\author{Yuqing Zhang}





