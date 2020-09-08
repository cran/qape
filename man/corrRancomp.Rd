\name{corrRancomp}
\alias{corrRancomp}
\title{Correction of predicted random components}


\description{The function computes the corrected predicted random components as presented in Chambers and Chandra (2013) in Section 2.2 to avoid the problem of underdispersion of residual bootstrap distributions.
}
\usage{
corrRancomp(model)
}
\arguments{
  \item{model}{ \emph{lmer}) object .
}
}
\value{the vector of corrected predicted random components.}

\references{
Chambers, R. and Chandra, H. (2013) A Random Effect Block Bootstrap for Clustered Data, Journal of Computational and Graphical Statistics, 22(2), 452-470.  \cr

}

\author{
Tomasz Zadlo}

\examples{

library(lme4)
data(invData)
attach(invData)
model=lmer(investments~newly_registered+(1|NUTS4))
corrRancomp(model)
}
