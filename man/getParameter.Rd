\name{getParameter}
\alias{getParameter}
\title{getParameter}
\description{Before the third step getting survival time, I need beta -- the coefficient 
from CoxBoost performed on the j=1,2,3... study. Here I think CoxBoost should 
be used on the original set instead of the one after resampling.
Also, grid, suvH, censH, which are useful for generating survival and censoring
times will be calculated here.
grid=grid corresponding to hazard estimations censH and survH
survH=cumulative Hazard for survival times distribution
censH=cumulative Hazard for censoring times distribution}
\usage{getParameter(esets, par_step)}
\arguments{
  \item{esets}{simulated ExpressionSets}
  \item{par_step}{CoxBoost parameter}
}

\value{returns a list of beta value for all the original data sets in order}

\author{Yuqing Zhang}





