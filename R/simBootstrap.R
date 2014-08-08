simBootstrap <- function
### the driver function to perform three-step bootstrap
(esets,
 ### the original ExpressionSets
 funSimData=simData,
 ### function to perform non-parametric bootstrap
 balance.variables=NULL,
 ### covariate names to balance in the simulated sets
 n.samples,
 ### number of samples to generate
 type="simulated",
 ### whether to include parametric bootstrap
 funParameter=getParameter,
 ### function to perform CoxBoost and return cumulative hazard of original sets
 step,
 ### step number of CoxBoost
 funSurvTime=generateSurvtime
 ### function to perform parametric bootstrap
){
  result <- funParameter(esets=esets, par_step=step)
  simmodels <- simData(esets=esets, balance.variables=balance.variables,
                       n.samples=n.samples, type=type)
  if(type=="simulated"){
    simmodels <- generateSurvtime(simmodels=simmodels, result=result)
  } 
  res <- list(esets.list=simmodels$esets, indices.list=simmodels$indices, setsID=simmodels$setsID, 
              lp.list=result$lp, beta.list=result$beta, 
              survH.list=result$survH, censH.list=result$censH, grid.list=result$grid)
  return(res)
  ### a list of values including:
  ### esets.list = a list of ExpressionSets, if type==simulated, then the survival time is
  ###              generated. 
  ### indices.list = a list of indices indicating which sample the simulated sample is in the 
  ###                original set
  ### setsID = a vector to indicate the original ID of simulated sets, if 
  ###          type=="original", setsID should be 1,2,3,...
  ### lp.list = a list of true linear predictor of each original data sets
  ### beta.list = a list of true coefficients used for simulating observations
}