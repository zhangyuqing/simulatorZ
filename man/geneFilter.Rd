\name{geneFilter}
\alias{geneFilter}
\title{geneFilter}
\description{the function to filter genes by Intergrative Correlation}
\usage{geneFilter(obj, cor.cutoff = 0.5)}
\arguments{
  \item{obj}{a list of ExpressionSet, matrix or SummarizedExperiment objects. If its elements are matrices,
columns represent samples, rows represent genes}
  \item{cor.cutoff}{the cutoff threshold for filtering genes. Only when the integrative correlation
between every pair of sets is larger than the cutoff value, will the gene 
be selected.}
}

\value{returns a list of ExpressionSets matrix or SummarizedExperiment objects
with genes filtered }
\references{Garrett-Mayer, E., Parmigiani, G., Zhong, X., Cope, L., 
Gabrielson, E., Cross-study validation and combined analysis of gene 
expression microarray data. Biostatistics. 2008 Apr;9(2):333-354.}
\author{Yuqing Zhang}




\examples{
set.seed(8)
library(curatedOvarianData)
library(GenomicRanges)
data(GSE17260_eset)
data(E.MTAB.386_eset)
data(GSE14764_eset)
esets <- list(GSE17260=GSE17260_eset, E.MTAB.386=E.MTAB.386_eset, GSE14764=GSE14764_eset)
esets.list <- lapply(esets, function(eset){
  return(eset[1:1500, 1:10])
})

result.set <- geneFilter(esets.list, 0)
result.set
### as we cannot calculate correlation with one set, this function just 
### delivers the same set if esets has length 1
result.oneset <- geneFilter(esets.list[1])
result.oneset

## Support matrices
X.list <- lapply(esets.list, function(eset){
  return(exprs(eset)) ## Columns represent samples!
})
result.set <- geneFilter(X.list, 0)
dim(result.set[[1]])

## Support SummarizedExperiment
nrows <- 200; ncols <- 6
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
rowRanges <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),
                     IRanges(floor(runif(200, 1e5, 1e6)), width=100),
                     strand=sample(c("+", "-"), 200, TRUE))
colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),
                     row.names=LETTERS[1:6])
sset <- SummarizedExperiment(assays=SimpleList(counts=counts),
                             rowRanges=rowRanges, colData=colData)
s.list <- list(sset, sset)
result.set <- geneFilter(s.list, 0.9) 
## the same set should resemble each other, no genes filtered
dim(assay(result.set[[1]]))
}
