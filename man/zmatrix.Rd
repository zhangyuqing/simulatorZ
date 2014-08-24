<<<<<<< HEAD
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





=======
\name{zmatrix}
\alias{zmatrix}
\title{zmatrix}
\description{generate a matrix of c statistics}
\usage{zmatrix(esets, y.vars, fold, trainingFun = masomenos, cvFun = funCV, 
    cvSubsetFun = cvSubsets)}
\arguments{
  \item{esets}{a list of ExpressionSets }
  \item{y.vars}{a list of response variables, all the response variables shold be
matrix, data.frame(with 2 columns) or Surv object}
  \item{fold}{cvFun parameter, in this case passes to funCV()}
  \item{trainingFun}{training function}
  \item{cvFun}{function to perform cross study within one set}
  \item{cvSubsetFun}{function to divide the expression sets into subsets for cross validation }
}

\value{outputs one matrix of validation statistics}

\author{Levi Waldron, Christoph Bernau, Yuqing Zhang}




\examples{
library(curatedOvarianData)
data( E.MTAB.386_eset )
eset1 <- E.MTAB.386_eset[1:100, 1:30]
eset2 <- E.MTAB.386_eset[1:100, 31:60]
eset3 <- E.MTAB.386_eset[1:100, 61:90]  
esets <- list(eset1, eset2, eset3) 

time1 <- eset1$days_to_death
#cens1 <- c(0, 0, 0, 1, 1)
cens1 <- sample(0:1, 30, replace=TRUE)
y1 <- Surv(time1, cens1)
time2 <- eset2$days_to_death
#cens2 <- c(1, 1, 0, 0, 0)
cens2 <- sample(0:1, 30, replace=TRUE)
y2 <- Surv(time2, cens2)
time3 <- eset3$days_to_death
#cens3 <- c(1, 0, 0, 0, 1)
cens3 <- sample(0:1, 30, replace=TRUE)
y3 <- Surv(time3, cens3)
y.vars <- list(y1, y2, y3)

# generate on original ExpressionSets
z <- zmatrix(esets, y.vars, 3)

# generate on simulated ExpressionSets
simmodels <- simBootstrap(esets, y.vars, 10, 100)
z <- zmatrix(simmodels$esets.list, simmodels$y.vars.list, 3)
z
}
>>>>>>> upstream/master
