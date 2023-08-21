\name{quantileNaN}
\alias{quantileNaN}

\title{quantile NaN}

\description{The function returns NaN when one of its arguments is NaN (instead of the error returned in this case by the classic quantile function)
 }

\usage{
quantileNaN(x, probs)
}

\arguments{
  \item{x}{numeric vector whose sample quantiles are wanted.}
  \item{probs}{numeric vector of probabilities with values in [0,1].}
}

\author{
Tomasz Zadlo}



