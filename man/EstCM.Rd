\name{EstCM}
\alias{EstCM}
\title{Estimated covariance matrix of predicted random effects
}
\description{A list of estimated covariance matrices of predicted random effects, where the length of the list equals the number of grouping variables used to define random effects as described in Carpenter, Goldstein and Rasbash (2003) in Section 3.2 and in Thai et al. (2013) in Section 2.3.3.
}
\usage{
EstCM(model)
}

\arguments{ \item{model}{ \emph{lmer}) object.
}
}
\value{list of estimated covariance matrices of predicted random effects.}

\references{
1. Carpenter, J.R., Goldstein, H. and Rasbash, J. (2003), A novel bootstrap procedure for assessing the relationship between class size and achievement. Journal of the Royal Statistical Society: Series C (Applied Statistics), 52, 431-443. \cr 

2. Thai, H.-T., Mentre, F., Holford, N.H., Veyrat-Follet, C. and Comets, E. (2013), A comparison of bootstrap approaches for estimating uncertainty of parameters in linear mixed-effects models. Pharmaceutical Statistics, 12, 129-140. \cr
}

\author{
Alicja Wolny-Dominiak, Tomasz Zadlo}

\examples{
library(lme4)
data(invData)
attach(invData)
model=lmer(investments~newly_registered+(1|NUTS4))
EstCM(model)
}

