\name{funTest}
\alias{funTest}
\title{funTest}
\description{test function}
\usage{funTest(eset, masomenosBeta, cstatFun)}
\arguments{
  \item{eset}{simulated ExpressionSet}
  \item{masomenosBeta}{return of masomenos(), optimal coefficient in Cox regression}
  \item{cstatFun}{function for calculating the c stat, in this case calcu_cstat()}
}

\value{returns the c statistics of cross study validation or cross validation}

\author{Yuqing Zhang}





