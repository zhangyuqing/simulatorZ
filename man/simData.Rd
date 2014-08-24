\name{simData}
\alias{simData}
\title{simData}
\description{simData() is a function to perform non-parametric bootstrap on a list of ExpressionSets. 
Input is a list of ExpressionSets(we will call it the original ones),  
output includes a list of the simulated ExpressionSets, as well as the indices of datasets and patients in the original sets. 
The function supports: 1.choosing whether to resample from sets by changing the parameter type, and 2. choosing whether to balance over part or all of the 
covariates by changing parameter balance.variables.}
\usage{simData(esets, n.samples, y.vars = list(), type = "two-steps", 
    balance.variables = NULL)}
\arguments{
  \item{esets}{a list of ExpressionSets, matrics or SummarizedExperiments. If elements are matricse,
columns represent samples}
  \item{n.samples}{an integer indicating how many samples should be resampled from each set}
  \item{y.vars}{a list of response variables}
  \item{type}{"one-step" or "two-steps". If type="one-step", the function will skip resampling the datasets, and directly resample from the original list
of ExpressionSets }
  \item{balance.variables}{balance.variables will be a vector of covariate names that should be balanced in the simulation. After balancing, the prevalence of covariate in each 
result set should be the same as the overall distribution across all original data sets. Default is set as NULL, when it will not balance over any covariate.
if isn't NULL, esets parameter should only be of class ExpressionSet}
}

\value{returns a list of simulated ExpressionSets, with names indicating its original set, and indices of the original patients.}

\author{Levi Waldron, Christoph Bernau, Yuqing Zhang}




\examples{
library(curatedOvarianData)
library(GenomicRanges)
data( E.MTAB.386_eset )
eset1 <- E.MTAB.386_eset[1:10, 1:5]
eset2 <- E.MTAB.386_eset[1:10, 6:10]
eset3 <- E.MTAB.386_eset[1:10, 11:15]
  
## simulate on multiple ExpressionSets
esets.list <- list(eset1, eset2, eset3)  
# one-step bootstrap: skip resampling set labels
simmodels <- simData(esets.list, 20, type="one-step")  
# two-step-non-parametric bootstrap
simmodels <- simData(esets.list, 10, type="two-steps")

## simulate one set
simmodels <- simData(list(eset1), 10, type="two-steps")

## balancing covariates
# single covariate
simmodels <- simData(list(eset2), 5, balance.variables="tumorstage")
# check the balancing effect
simmodels$prob.desired
table(eset2$tumorstage) / sum(table(eset2$tumorstage))
table(simmodels$esets[[1]]$tumorstage) / sum( table(simmodels$esets[[1]]$tumorstage))
# multiple covariates
simmodels <- simData(list(eset2), 5, 
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
