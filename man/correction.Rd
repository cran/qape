\name{correction}
\alias{correction}
\title{Correction term for predicted random effects
}
\description{The function computes the list of matrices used to correct predicted random effects as presented in Carpenter, Goldstein and Rasbash (2003) in Section 3.2 to avoid the problem of underdispersion of residual bootstrap distributions.
}
\usage{
correction(model)
}

\arguments{ \item{model}{ \emph{lmer} object.
}
}

\value{a list of square matrices used to correct predicted random effects. The length of the list is equal the number of grouping variables used in case of random effects. Each matrix is of order equal the number of random effects at the considered level of grouping.
}

\references{
Carpenter, J.R., Goldstein, H. and Rasbash, J. (2003), A novel bootstrap procedure for assessing the relationship between class size and achievement. Journal of the Royal Statistical Society: Series C (Applied Statistics), 52, 431-443. \cr
}

\author{
Tomasz Zadlo}

\examples{
library(lme4)
data(invData)
attach(invData)
model <- lmer(investments ~ newly_registered + ((1|NUTS2) + 
((newly_registered - 1)|NUTS2) + ((newly_registered)|NUTS4)))
correction(model)

detach(invData)
}
