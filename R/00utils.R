### a generic function with methods for matrix, ExpressionSet, and RangedSummarizedExperiment objects.
### Matrix, ExpressionSet or RangedSummarizedExperiment

getMatrix <- function(obj){
    if (is.matrix(obj))
        return(obj)
    if (is(obj, "ExpressionSet"))
        return(exprs(obj))
    if (is(obj, "RangedSummarizedExperiment")
     || is(obj, "SummarizedExperiment"))  # for backward compatibility
        return(assay(obj))
	### return the expression matrix
    stop("Wrong class of obj!")
}


changevar <- function(dat, var_names){
  ### change discrete/categorical variables into binary variables
  if(!is.data.frame(dat)){dat <- as.data.frame(dat)}
  colnames(dat) <- var_names
  
  new_var_names <- c()
  res_dat <- matrix(0, nrow=nrow(dat), ncol=0, dimnames=list(rownames(dat), NULL))  
  
  for(i in seq_along(colnames(dat))){
    if(class(dat[, i])=="factor"){tmp <- levels(dat[, i])}else{tmp <- unique(dat[, i])}
    tmp <- tmp[order(tmp)]
    
    if(length(tmp)<=2){
      res_dat <- cbind(res_dat, model.matrix(~as.factor(dat[,i]))[,2])
      new_var_names <- c(new_var_names, var_names[i])
    }else{
      tmp_dat <- matrix(0, nrow=nrow(dat), ncol=length(tmp)-1)
      for(j in 1:(length(tmp)-1)){
        ii <- which(dat[,i]==tmp[j])
        tmp_dat[ii, j] <- 1
        tmp_dat[-ii, j] <- 0
      }
      
      new_var_names <- c(new_var_names, paste(colnames(dat)[i], 1:(length(tmp)-1), sep="_"))
      res_dat <- cbind(res_dat, tmp_dat)
    }
  }
  
  colnames(res_dat) <- new_var_names
  return(res_dat)
  ### a data frame with all variables binarized
}

