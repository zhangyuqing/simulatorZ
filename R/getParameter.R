getParameter <- function
### Before the third step getting survival time, I need beta -- the coefficient 
### from CoxBoost performed on the j=1,2,3... study. Here I think CoxBoost should 
### be used on the original set instead of the one after resampling.
### Also, grid, suvH, censH, which are useful for generating survival and censoring
### times will be calculated here.
### grid=grid corresponding to hazard estimations censH and survH
### survH=cumulative Hazard for survival times distribution
### censH=cumulative Hazard for censoring times distribution
(esets,
 ### simulated ExpressionSets
 par_step 
 ### CoxBoost parameter
 #par_penalty
 ### CoxBoost parameter
){
  Beta <- grid <- survH <- censH <- result <- lp <- list()
  for(i in 1:length(esets)){
    print(i)
    ### PART 1: get beta
    time <- esets[[i]]$dmfs.time
    time <- as.numeric(as.character(time))
    status <- esets[[i]]$dmfs.cens
    status <- as.numeric(as.character(status))
    cbfit <- CoxBoost(time=time, status=status, x=t(exprs(esets[[i]])), 
                      stepno=par_step, standardize=FALSE)
    Beta[[i]] <- coef(cbfit)    
    ### PART 2: get grid, suvH, censH
    ### Survival time ==> dmfs.cens=1, Censoring time ==> dmfs.cens=0
    ### get survH, censH, grid
    X <- exprs(esets[[i]])
    lp[[i]] <- as.numeric(Beta[[i]] %*% X)  ## Calculate linear predictor
    grid[[i]] <- seq(0, max(time), by = 1)
    survH[[i]] <- basehaz.gbm(t=time, delta=status, f.x=lp[[i]], 
                              t.eval=grid[[i]], smooth=TRUE, cumulative=TRUE)    
    inverse_status <- (-1) * status + 1
    censH[[i]] <- basehaz.gbm(t=time, delta=inverse_status, f.x=rep(0, length(time)), 
                              t.eval=grid[[i]], smooth=TRUE, cumulative=TRUE)
  }
  result <- list(beta=Beta, grid=grid, survH=survH, censH=censH, lp=lp)
  return(result)
  ### returns a list of beta value for all the original data sets in order
}