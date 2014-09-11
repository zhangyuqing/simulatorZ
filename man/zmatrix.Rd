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
\usage{zmatrix(obj, y.vars, fold, trainingFun = masomenos, cvFun = funCV, 
    cvSubsetFun = cvSubsets)}
\arguments{
  \item{obj}{a list of ExpressionSet, matrix or SummarizedExperiment objects, if its 
elements are matrices, columns represent samples }
  \item{y.vars}{a list of response variables, all the response variables shold be
matrix, data.frame(with 2 columns) or Surv object}
  \item{fold}{cvFun parameter, in this case passes to funCV()}
  \item{trainingFun}{training function}
  \item{cvFun}{function to perform cross study within one set}
  \item{cvSubsetFun}{function to divide the expression sets into subsets for cross validation }
}

\value{outputs one matrix of validation statistics}

\author{Yuqing Zhang, Christoph Bernau, Levi Waldron}




\examples{
library(curatedOvarianData)
library(GenomicRanges)
source(system.file("extdata", "patientselection.config",
                   package="curatedOvarianData"))
source(system.file("extdata", "createEsetList.R", package="curatedOvarianData"))
esets.list <- lapply(esets, function(eset){
  return(eset[1:500, 1:20])
})
esets.list <- esets.list[1:5]
## simulate on multiple ExpressionSets
set.seed(8) 

y.list <- lapply(esets.list, function(eset){
  time <- eset$days_to_death
  cens.chr <- eset$vital_status
  cens <- c()
  for(i in 1:length(cens.chr)){
    if(cens.chr[i] == "living") cens[i] <- 1
    else cens[i] <- 0
  }
  y <- Surv(time, cens)
  return(y)
})

# generate on original ExpressionSets
z <- zmatrix(esets.list, y.list, 3)

# generate on simulated ExpressionSets
simmodels <- simBootstrap(esets.list, y.list, 10, 100)
z <- zmatrix(simmodels$obj.list, simmodels$y.vars.list, 3)

# support matrix
X.list <- lapply(esets.list, function(eset){
  newx <- exprs(eset) ### columns represent samples !!
  return(newx)
}) 
z <- zmatrix(X.list, y.list, 3)

# support SummarizedExperiment
nrows <- 200; ncols <- 6
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
rowData <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),
                   IRanges(floor(runif(200, 1e5, 1e6)), width=100),
                   strand=sample(c("+", "-"), 200, TRUE))
colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),
                     row.names=LETTERS[1:6])
sset <- SummarizedExperiment(assays=SimpleList(counts=counts),
                             rowData=rowData, colData=colData)

time <- sample(4500:4700, 6, replace=TRUE)
cens <- sample(0:1, 6, replace=TRUE)
y.vars <- Surv(time, cens)

z <- zmatrix(list(sset[,1:3], sset[,4:6]), list(y.vars[1:3,],y.vars[4:6,]), 3)
}
>>>>>>> upstream/master
