\name{Zfun}
\alias{Zfun}

\title{Matrix Z creator}

\description{The function creates the Z matrix of auxiliary variables 
asscociatied with random effects. 
}
\usage{
Zfun(model, data)
}

\arguments{
  \item{model}{formula of model (use \emph{formula()} function).
}
  \item{data}{data.
}
}
\value{
\item{Z}{Z matrix.}
\item{vNames}{labels of random effects.}
\item{ZBlockNames}{labels of blocks of random effects.}
}

\author{
Alicja Wolny-Dominiak}

\examples{
data(invData)
modelFormula <- formula(investments~newly_registered + (newly_registered | NUTS2))
reg <- invData

Zfun(modelFormula, reg)
} 


