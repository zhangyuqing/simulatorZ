simBootstrap <- structure(function
### the driver function to perform three-step bootstrap resampling 
### to get independent genomic data sets
(obj,
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
  result <- funTrueModel(obj=obj, y.vars=y.vars, parstep=parstep)
  simmodels <- funSimData(obj=obj, balance.variables=balance.variables,
                       n.samples=n.samples, type=type, y.vars=y.vars)
  simmodels <- funSurvTime(simmodels=simmodels, result=result) 
  res <- list(obj.list=simmodels$obj, y.vars.list=simmodels$y.vars,
              indices.list=simmodels$indices, setsID=simmodels$setsID, 
              lp.list=result$lp, beta.list=result$beta, 
              survH.list=result$survH, censH.list=result$censH, grid.list=result$grid)
  return(res)
  ### a list of values including:
  ### obj.list = a list of simulated objects the same type as input
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
    for(i in 1:length(cens.chr)){
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
  simmodels <- simBootstrap(obj=s.list, y.vars=y.vars, 100, 100) 
})