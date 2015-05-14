\name{simBootstrap}
\alias{simBootstrap}
\title{simBootstrap}
\description{the driver function to perform three-step bootstrap resampling 
to get independent genomic data sets}
\usage{simBootstrap(obj, y.vars, n.samples, parstep, type = "two-steps", 
    balance.variables = NULL, funSimData = simData, funTrueModel = getTrueModel, 
    funSurvTime = simTime)}
\arguments{
  \item{obj}{a list of ExpressionSet, matrix or RangedSummarizedExperiment  }
  \item{y.vars}{a list of reponse variables, elements can be class Surv, matrix or data.frame}
  \item{n.samples}{number of samples to resample in each set}
  \item{parstep}{step number to fit CoxBoost}
  \item{type}{whether to include resampling set labels }
  \item{balance.variables}{covariate names to balance in the simulated sets}
  \item{funSimData}{function to perform non-parametric bootstrap}
  \item{funTrueModel}{function to construct true models in original sets}
  \item{funSurvTime}{function to perform parametric bootstrap}
}

\value{a list of values including:
obj.list = a list of simulated objects the same type as input
indices.list = a list of indices indicating which sample the simulated sample is in the 
original set
setsID = a vector to indicate the original ID of simulated sets, if 
type=="original", setsID should be 1,2,3,...
lp.list = a list of true linear predictor of each original data sets
beta.list = a list of true coefficients used for simulating observations
survH.list = list of cumulative survival hazard
censH.list = list of cumulative censoring hazard
grid.list = list of timeline grid corresponding to survH and censH respectivley}

\author{Yuqing Zhang, Levi Waldron}




\examples{
library(curatedOvarianData)
library(GenomicRanges)
data(GSE17260_eset)
data(E.MTAB.386_eset)
data(GSE14764_eset)
esets <- list(GSE17260=GSE17260_eset, E.MTAB.386=E.MTAB.386_eset, GSE14764=GSE14764_eset)
esets.list <- lapply(esets, function(eset){
  return(eset[1:500, 1:20])
})

## simulate on multiple ExpressionSets
set.seed(8) 

y.list <- lapply(esets.list, function(eset){
  time <- eset$days_to_death
  cens.chr <- eset$vital_status
  cens <- c()
  for(i in seq_along(cens.chr)){
    if(cens.chr[i] == "living") cens[i] <- 1
    else cens[i] <- 0
  }
  y <- Surv(time, cens)
  return(y)
})

simmodels <- simBootstrap(obj=esets.list, y.vars=y.list, 10, 100)
simmodels$obj.list[[1]]

# balance covariates
simmodels <- simBootstrap(obj=esets.list, y.vars=y.list, 10, 100,
                          balance.variables="tumorstage")

## Support RangedSummarizedExperiment
nrows <- 200; ncols <- 10
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
rowRanges <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),
                     IRanges(floor(runif(200, 1e5, 1e6)), width=100),
                     strand=sample(c("+", "-"), 200, TRUE))
colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 5),
                     row.names=LETTERS[1:10])
sset <- SummarizedExperiment(assays=SimpleList(counts=counts),
                             rowRanges=rowRanges, colData=colData)

s.list <- list(sset[,1:5], sset[,6:10])
time <- c(540, 527, 668, 587, 620, 540, 527, 668, 587, 620)
cens <- c(1, 0, 0, 1, 0, 1, 0, 0, 1, 0)
y.vars <- Surv(time, cens)
y.vars <- list(y.vars[1:5,],y.vars[1:5,])
simmodels <- simBootstrap(obj=s.list, y.vars=y.vars, 100, 100) 
}
