masomenos <- structure(function
### function for Mas-o-menos algorithm
(X, 
 ### gene expression matrix, with rows corresponding to patients and columns corresponding to genes 
 y 
 ### response variable, a data.frame, matrix or Surv object, with the first column as time to event and the second column 
 ### censoring status 
){
  n.genes <- ncol(X)
  beta <- numeric(n.genes) 
  X <- data.frame(X, row.names=NULL)  
  y <- Surv(y[, 1], y[, 2])
  alpha <- rowCoxTests(t(X), y, option="fast")
  for(p in 1:(n.genes)){
    ## Calculate v
    if(alpha[p, 1] > 0) beta[p] <- 1
    if(alpha[p, 1] == 0) beta[p] <- 0
    if(alpha[p, 1] < 0) beta[p] <- -1
    beta[p] <- beta[p] / n.genes
  }
  return(beta)
  ### return betas
},ex=function(){
  library(curatedOvarianData)
  data( E.MTAB.386_eset )
  eset <- E.MTAB.386_eset[1:100, 1:30]
  
  X <- t(exprs(eset))
  
  time <- eset$days_to_death
  cens <- sample(0:1, 30, replace=TRUE)
  y <- Surv(time, cens)
  
  beta <- masomenos(X, y)
  beta
})