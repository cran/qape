\name{corrRanef}
\alias{corrRanef}
\title{Correction of predicted random effects}

\description{The function computes the corrected predicted random effects as presented in Carpenter, Goldstein and Rasbash (2003) in Section 3.2 to avoid the problem of underdispersion of residual bootstrap distributions.
}
\usage{
corrRanef(model)
}
\arguments{
\item{model}{ \emph{lmer}) object .
}
}
\value{a list of corrected predicted random effects (of the same form as \emph{ranef(model)})}

\references{
Carpenter, J.R., Goldstein, H. and Rasbash, J. (2003), A novel bootstrap procedure for assessing the relationship between class size and achievement. Journal of the Royal Statistical Society: Series C (Applied Statistics), 52, 431-443. \cr
}

\author{
Tomasz Zadlo}

\examples{
library(lme4)
data(invData)
attach(invData)
model=lmer(investments~newly_registered+(1|NUTS4))
corrRanef(model)
}

