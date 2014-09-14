\name{rowCoxTests}
\alias{rowCoxTests}
\title{rowCoxTests}
\description{method for performing Cox regression}
\usage{rowCoxTests(X, y, option = c("fast", "slow"), ...)}
\arguments{
  \item{X}{Gene expression data. The following formats are available:   
matrix Rows correspond to observations, columns to variables.
data.frame Rows correspond to observations, columns to variables.
ExpressionSet rowCoxTests will extract the expressions using exprs().}
  \item{y}{Survival Response, an object of class:
Surv if X is of type data.frame or matrix
character if X is of type ExpressionSet. 
In this case y is the name of the survival 
response in the phenoData of X. If survival 
time and indicator are stored separately 
in the phenoData one can specify a two-element 
character vector the first element representing 
the survival time variable.}
  \item{option}{"fast" loops over rows in C, "slow" calls coxph 
directly in R. The latter method may be used if 
something goes wrong with the "fast" method.}
  \item{\dots}{currently unused}
}

\value{dataframe with two columns: coef = Cox regression
coefficients, p.value =
Wald Test p-values.  Rows correspond to the rows of X.}

\author{Christoph Bernau}




\examples{
#test
##regressor-matrix (gene expressions)
X<-matrix(rnorm(1e6),nrow=10000)
#seed
set.seed(123)
#times
time<-rnorm(n=ncol(X),mean=100)
#censoring(1->death)
status<-rbinom(n=ncol(X),size=1, prob=0.8)

##survival object
y<-Surv(time,status)

## Do 10,000 Cox regressions:
system.time(output <- rowCoxTests(X=X,y=y, option="fast"))
}
