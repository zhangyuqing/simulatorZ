masomenos <- function
### Mas-o-menos algorithm, mainly used as training function
(X, 
 ### mostly assayData of an ExpressionSet object 
 y, 
 ### response, using two columns "dmfs.time" and "dmfs.cens" in the phenoData of an ExpressionSet object
 n.genes,
 ### number of gene expression variables
 p.cutoff 
 ### p threshold ??
){ 
  print("TRAINING")
  beta <- numeric(n.genes) 
  X <- data.frame(X, row.names=NULL)  
  y <- Surv(y[, 1], y[, 2])
  alpha <- rowCoxTests(t(X), y, option="fast")
  for(p in 1:(n.genes)){
    ## Calculate v
    if(alpha[p, 1] > 0) beta[p] <- 1
    if(alpha[p, 1] == 0) beta[p] <- 0
    if(alpha[p, 1] < 0) beta[p] <- -1
    beta[p] <- beta[p] / sqrt(n.genes)
  }
  return(alpha[, 1])
  ### return betas
}