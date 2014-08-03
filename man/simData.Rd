\name{simData}
\alias{simData}
\title{simData}
\description{Input is a list of ExpressionSets of the original sets, 
output is another list of ExpressionSets of the simulated sets
Contains the two steps of non-parametric bootstrap }
\usage{simData(esets, balance.variables)}
\arguments{
  \item{esets}{original data sets}
  \item{balance.variables}{Later, balance.variables will be a vector of variable names that should be balanced in the simulation.}
}

\value{returns a list of simulated ExpressionSets, with names indicating its origin.}

\author{Yuqing Zhang}





