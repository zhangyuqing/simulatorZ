\name{funCV}
\alias{funCV}
\title{funCV}
\description{Cross validation function}
\usage{funCV(eset, fold, p.cutoff, trainFun, testFun, cstatFun)}
\arguments{
  \item{eset}{simulated ExpressionSets }
  \item{fold}{tells the number of folds in cross validation}
  \item{p.cutoff}{trainFun parameter, in this case passes to masomenos()}
  \item{trainFun}{training function, in this case masomenos()}
  \item{testFun}{test function, in this case funTest()}
  \item{cstatFun}{function for calculating the c stat, in this case calcu_cstat()}
}

\value{returns the c statistics of cross validation(CV)}

\author{Yuqing Zhang}





