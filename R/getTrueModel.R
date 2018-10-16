getTrueModel <- structure(function
### The parametric bootstrap simulation depends on the true model of original sets.
### This function is to generate useful values from the true models for further analysis.  
### We fit CoxBoost to the original sets and use the coefficients to simulate
### the survival and censoring time. grid, survH, censH, which are useful for this purpose.
### grid=grid corresponding to hazard estimations censH and survH
### survH=cumulative hazard for survival times distribution
### censH=cumulative hazard for censoring times distribution
(esets,
 ### a list of ExpressionSets, matrix or SummarizedExperiment
 y.vars,
 ### a list of response variables
 parstep, 
 ### CoxBoost parameter
 #par_penalty
 ### CoxBoost parameter
 balance.variables=NULL
 ### variable names to be balanced.
){
  Beta <- grid <- survH <- censH <- result <- lp <- list()
  for(i in seq_along(esets)){
    print(i)
    ### PART 1: get beta
    time <- y.vars[[i]][, 1]
    time <- as.numeric(as.character(time))
    status <- y.vars[[i]][, 2]
    status <- as.numeric(as.character(status))
    
    X <- t(getMatrix(esets[[i]])) 
    
    unpen.ind <- NULL  
    if(!is.null(balance.variables)){
      X1 <- pData(esets[[i]])[,balance.variables]
      X2 <- changevar(X1, balance.variables)
      X <- cbind(X2, X)
      unpen.ind <- 1:ncol(X2)
      
      if(ncol(X2)==1){
        varseq <- var(X2[which(status==1),])
        if(varseq==0){unpen.ind <- NULL}
      }else{
        varseq <- apply(X2[which(status==1),],2,var)
        unpen.ind <- unpen.ind[which(varseq!=0)]
      }
    }
    
    cbfit <- CoxBoost(time=time, status=status, x=X, unpen.index=unpen.ind,
                      stepno=parstep, standardize=FALSE)
    
    Beta[[i]] <- coef(cbfit)    
    
    ### PART 2: get grid, suvH, censH
    ### Survival time ==> dmfs.cens=1, Censoring time ==> dmfs.cens=0
    ### get survH, censH, grid    
    
    lp[[i]] <- as.numeric(Beta[[i]] %*% t(X))  ## Calculate linear predictor
    grid[[i]] <- seq(0, max(time), by = 1)
    
    survH[[i]] <- basehaz.gbm(t=time, delta=status, f.x=lp[[i]], 
                              t.eval=grid[[i]], smooth=TRUE, cumulative=TRUE) 
    
    inverse_status <- (-1) * status + 1
    censH[[i]] <- basehaz.gbm(t=time, delta=inverse_status, f.x=rep(0, length(time)), 
                              t.eval=grid[[i]], smooth=TRUE, cumulative=TRUE)
    
  }
  result <- list(beta=Beta, grid=grid, survH=survH, censH=censH, lp=lp)  
  return(result)
  ### returns a list of values:
  ### beta: True coefficients obtained by fitting CoxBoost to the original ExpressionSets
  ### grid: timeline grid corresponding to hazard estimations censH and survH
  ### survH: cumulative hazard for survival times distribution
  ### censH: cumulative hazard for censoring times distribution
  ### lp: true linear predictors 
},ex=function(){
  library(curatedOvarianData)
  data(GSE14764_eset)
  data(E.MTAB.386_eset)
  esets.list <- list(GSE14764=GSE14764_eset[1:500, 1:20], 
                     E.MTAB.386=E.MTAB.386_eset[1:500, 1:20])
  rm(E.MTAB.386_eset, GSE14764_eset)
  
  ## simulate on multiple ExpressionSets
  set.seed(8) 
  
  y.list <- lapply(esets.list, function(eset){
    time <- eset$days_to_death
    cens.chr <- eset$vital_status
    cens <- rep(0, length(cens.chr))
    cens[cens.chr=="living"] <- 1
    return(Surv(time, cens))
  })
     
  res1 <- getTrueModel(esets.list, y.list, 100)
  ## Get true model from one set
  res2 <- getTrueModel(esets.list[1], y.list[1], 100)
  names(res2)
  res2$lp
  ## note that y.list[1] cannot be replaced by y.list[[1]]
})
