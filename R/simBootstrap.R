simBootstrap <- structure(function
### the driver function to perform three-step bootstrap simulation
(esets,
 ### a list of ExpressionSet, matrix or SummarizedExperiment  
 y.vars,
 ### a list of reponse variables, elements can be class Surv, matrix or data.frame
 n.samples,
 ### number of samples to resample in each set
 parstep,
 ### step number to fit CoxBoost
 type="two-steps",
 ### whether to include resampling set labels 
 balance.variables=NULL,
 ### covariate names to balance in the simulated sets
 funSimData=simData,
 ### function to perform non-parametric bootstrap
 funTrueModel=getTrueModel,
 ### function to construct true models in original sets
 funSurvTime=simTime
 ### function to perform parametric bootstrap
){
  result <- funTrueModel(esets=esets, y.vars=y.vars, parstep=parstep)
  simmodels <- funSimData(esets=esets, balance.variables=balance.variables,
                       n.samples=n.samples, type=type, y.vars=y.vars)
  simmodels <- funSurvTime(simmodels=simmodels, result=result) 
  res <- list(esets.list=simmodels$esets, y.vars.list=simmodels$y.vars,
              indices.list=simmodels$indices, setsID=simmodels$setsID, 
              lp.list=result$lp, beta.list=result$beta, 
              survH.list=result$survH, censH.list=result$censH, grid.list=result$grid)
  return(res)
  ### a list of values including:
  ### esets.list = a list of simulated ExpressionSets
  ### indices.list = a list of indices indicating which sample the simulated sample is in the 
  ###                original set
  ### setsID = a vector to indicate the original ID of simulated sets, if 
  ###          type=="original", setsID should be 1,2,3,...
  ### lp.list = a list of true linear predictor of each original data sets
  ### beta.list = a list of true coefficients used for simulating observations
  ### survH.list = list of cumulative survival hazard
  ### censH.list = list of cumulative censoring hazard
  ### grid.list = list of timeline grid corresponding to survH and censH respectivley
  
},ex=function(){
  library(curatedOvarianData)
  library(GenomicRanges)
  data( E.MTAB.386_eset )
  eset1 <- E.MTAB.386_eset[1:10, 1:5]
  eset2 <- E.MTAB.386_eset[1:10, 6:10]
  eset3 <- E.MTAB.386_eset[1:10, 11:15]  
  esets <- list(eset1, eset2, eset3) 
  
  time1 <- eset1$days_to_death
  cens1 <- c(0, 0, 0, 1, 1)
  y1 <- Surv(time1, cens1)
  time2 <- eset2$days_to_death
  cens2 <- c(1, 1, 0, 0, 0)
  y2 <- Surv(time2, cens2)
  time3 <- eset3$days_to_death
  cens3 <- c(1, 0, 0, 0, 1)
  y3 <- Surv(time3, cens3)
  y.vars <- list(y1, y2, y3)
  
  simmodels <- simBootstrap(esets=esets, y.vars=y.vars, 10, 100)
  simmodels$esets.list[[1]]
  
  # skip resampling labels
  simmodels <- simBootstrap(esets=esets, y.vars=y.vars, 10, 100,
                            type="one-step")
  
  # balance covariates
  simmodels <- simBootstrap(esets=esets, y.vars=y.vars, 10, 100,
                            balance.variables="tumorstage")
  
  ## Support matrices
  X.list <- lapply(esets, function(eset){
    return(exprs(eset))
  })
  simmodels <- simBootstrap(esets=X.list, y.vars=y.vars, 10, 100)
  
  ## Support SummarizedExperiment
  nrows <- 200; ncols <- 10
  counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
  rowData <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),
                     IRanges(floor(runif(200, 1e5, 1e6)), width=100),
                     strand=sample(c("+", "-"), 200, TRUE))
  colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 5),
                       row.names=LETTERS[1:10])
  sset <- SummarizedExperiment(assays=SimpleList(counts=counts),
                               rowData=rowData, colData=colData)
  
  s.list <- list(sset[,1:5], sset[,6:10])
  time <- c(540, 527, 668, 587, 620, 540, 527, 668, 587, 620)
  cens <- c(1, 0, 0, 1, 0, 1, 0, 0, 1, 0)
  y.vars <- Surv(time, cens)
  y.vars <- list(y.vars[1:5,],y.vars[1:5,])
  simmodels <- simBootstrap(esets=s.list, y.vars=y.vars, 100, 100) 
})