<<<<<<< HEAD
\name{masomenos}
\alias{masomenos}
\title{masomenos}
\description{Mas-o-menos algorithm, mainly used as training function}
\usage{masomenos(X, y, n.genes, p.cutoff)}
\arguments{
  \item{X}{mostly assayData of an ExpressionSet object }
  \item{y}{response, using two columns "dmfs.time" and "dmfs.cens" in the phenoData of an ExpressionSet object}
  \item{n.genes}{number of gene expression variables}
  \item{p.cutoff}{p threshold ??}
}

\value{return betas}

\author{Yuqing Zhang}





=======
\name{masomenos}
\alias{masomenos}
\title{masomenos}
\description{function for Mas-o-menos algorithm}
\usage{masomenos(X, y)}
\arguments{
  \item{X}{gene expression matrix, with rows corresponding to patients and columns corresponding to genes }
  \item{y}{response variable, a data.frame, matrix or Surv object, with the first column as time to event and the second column 
censoring status }
}

\value{return betas}

\author{Levi Waldron, Christoph Bernau, Yuqing Zhang}




\examples{
library(curatedOvarianData)
data( E.MTAB.386_eset )
eset <- E.MTAB.386_eset[1:100, 1:30]

X <- t(exprs(eset))

time <- eset$days_to_death
cens <- sample(0:1, 30, replace=TRUE)
y <- Surv(time, cens)

beta <- masomenos(X, y)
beta
}
>>>>>>> upstream/master
