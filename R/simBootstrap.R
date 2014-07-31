simBootstrap <- function
### the driver function to perform three-step bootstrap
(esets,
 ### the original ExpressionSets
 funSimData,
 ### function to perform non-parametric bootstrap
 balance.variables,
 ### covariate names to balance in the simulated sets
 funParameter,
 ### function to perform CoxBoost and return cumulative hazard of original sets
 par_step,
 ### step number of CoxBoost
 funSurvTime
 ### function to perform parametric bootstrap
){
  result <- funParameter(esets=esets, par_step=par_step)
  new_esets <- simData(esets=esets, balance.variables=balance.variables)  
  new_esets <- generateSurvtime(esets=new_esets, result=result)
  return(new_esets)
  ### a list of simulated ExpressionSets
}