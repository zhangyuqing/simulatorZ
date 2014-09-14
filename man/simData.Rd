\name{simData}
\alias{simData}
\title{simData}
\description{simData is a function to perform non-parametric bootstrap resampling
on a list of (original) data sets, both on set level and patient level,
in order to simulate independent genomic sets. }
\usage{simData(obj, n.samples, y.vars = list(), type = "two-steps", 
    balance.variables = NULL)}
\arguments{
  \item{obj}{a list of ExpressionSets, matrics or SummarizedExperiments. If elements are 
matricse, columns represent samples}
  \item{n.samples}{an integer indicating how many samples should be resampled from each set}
  \item{y.vars}{a list of response variables, can be Surv object, or matrix or data.frame
with two columns}
  \item{type}{string "one-step" or "two-steps". If type="one-step", the function will 
skip resampling the datasets, and directly resample from the original list
of obj }
  \item{balance.variables}{balance.variables will be a vector of covariate names that should be 
balanced in the simulation. After balancing, the prevalence of covariate 
in each result set should be the same as the overall distribution across 
all original data sets. Default is set as NULL, when it will not balance 
over any covariate. if isn't NULL, esets parameter should only be of class
ExpressionSet}
}

\value{returns a list of simulated ExpressionSets, with names indicating its original set, and indices of the original patients.
prob.desired and prob.real are only useful when balance.varaibles is set.
prob.desired shows overall distrubition of the specified covariate. prob.list
shows the sampling probability in each set after balancing}

\author{Yuqing Zhang, Levi Waldron}




\examples{
library(curatedOvarianData)
library(GenomicRanges)

data(GSE17260_eset)
data(E.MTAB.386_eset)
data(GSE14764_eset)

esets <- list(GSE17260=GSE17260_eset, E.MTAB.386=E.MTAB.386_eset, GSE14764=GSE14764_eset)
esets.list <- lapply(esets, function(eset){
  return(eset[1:1000, 1:10])
})

## simulate on multiple ExpressionSets
set.seed(8)
# one-step bootstrap: skip resampling set labels
simmodels <- simData(esets.list, 20, type="one-step")  
# two-step-non-parametric bootstrap
simmodels <- simData(esets.list, 10, type="two-steps")

## simulate one set
simmodels <- simData(list(esets.list[[1]]), 10, type="two-steps")

## balancing covariates
# single covariate
simmodels <- simData(list(esets.list[[1]]), 5, balance.variables="tumorstage")

# multiple covariates
simmodels <- simData(list(esets.list[[1]]), 5, 
                     balance.variables=c("tumorstage", "age_at_initial_pathologic_diagnosis"))  

## Support matrices
X.list <- lapply(esets.list, function(eset){
  return(exprs(eset))
})
simmodels <- simData(X.list, 20, type="two-steps")

## Support SummarizedExperiment
nrows <- 200; ncols <- 6
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
rowData <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),
                   IRanges(floor(runif(200, 1e5, 1e6)), width=100),
                   strand=sample(c("+", "-"), 200, TRUE))
colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),
                     row.names=LETTERS[1:6])
sset <- SummarizedExperiment(assays=SimpleList(counts=counts),
                             rowData=rowData, colData=colData)

s.list <- list(sset[,1:3], sset[,4:6])
simmodels <- simData(s.list, 20, type="two-steps")
}
