funTest <- function
### test function
(eset,
 ### simulated ExpressionSet
 masomenosBeta, 
 ### return of masomenos(), optimal coefficient in Cox regression
 cstatFun
 ### function for calculating the c stat, in this case calcu_cstat()
){
  print("TESTING")
  Xj <- t(exprs(eset))
  Xj <- scale(Xj)
  Sj <- Xj %*% masomenosBeta
  z <- cstatFun(eset, Sj)
  return(z)
  ### returns the c statistics of cross study validation or cross validation
}