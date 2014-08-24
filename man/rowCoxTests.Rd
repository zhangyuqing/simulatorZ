\name{rowCoxTests}
\alias{rowCoxTests}
\title{rowCoxTests}
\description{method for performing Cox regression}
\usage{rowCoxTests(X, y, option = c("fast", "slow"), ...)}
\arguments{
  \item{X}{
}
  \item{y}{
}
  \item{option}{
}
  \item{\dots}{
}
}

\value{dataframe with two columns: coef = Cox regression
coefficients, p.value =
Wald Test p-values.  Rows correspond to the rows of X.}

\author{Levi Waldron, Christoph Bernau, Yuqing Zhang}





