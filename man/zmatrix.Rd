\name{zmatrix}
\alias{zmatrix}
\title{zmatrix}
\description{Takes a list of ExpressionSets to generate one matrix of validation statistics
already-defined training, test and CV functions. }
\usage{zmatrix(esets, y.vars, trainingFun, testFun, cvFun, cstatFun, 
    p.cutoff, fold)}
\arguments{
  \item{esets}{simulated ExpressionSets }
  \item{y.vars}{some strings to indicate the response variable}
  \item{trainingFun}{training function, in this case masomenos()}
  \item{testFun}{test function, in this case funTest()}
  \item{cvFun}{function to do cross validation, in this case funCV()}
  \item{cstatFun}{function to calculate the c statistics, in this case calcu_cstat}
  \item{p.cutoff}{trainFun parameter, in this case passes to masomenos()}
  \item{fold}{cvFun parameter, in this case passes to funCV()}
}

\value{outputs one matrix of validation statistics}

\author{Yuqing Zhang}





