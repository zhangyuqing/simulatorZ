\name{generateSurvtime}
\alias{generateSurvtime}
\title{generateSurvtime}
\description{Step.3 Generating survival time --- parametric bootstrap
update the dmfs.time in simulated ExpressionSets}
\usage{generateSurvtime(esets, result)}
\arguments{
  \item{esets}{simulated ExpressionSets }
  \item{result}{return of getParameter(), a list which consists of Beta, grid, survH and censH
for each original data set}
}

\value{survival time is saved in phenodata, here the function still returns the ExpressionSets}

\author{Yuqing Zhang}





