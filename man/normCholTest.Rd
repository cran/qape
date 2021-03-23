\name{normCholTest}
\alias{normCholTest}

\title{Test of normality of the dependent variable}

\description{The function conducts a test of normality of the depenedent 
variable based on residuals transformed using Cholesky decomposition of
the inverse of the estimated variance-covariance matrix of the variable.
}
\usage{
normCholTest(model, normTest)
}

\arguments{
  \item{model}{lmer object.
}
  \item{normTest}{function which implements a normality test e.g. shapiro.test
  (takes a vector of the values of the variable as an argument and 
  conducts a test of normality of the variable).
}
}
\value{
\item{testResults}{output of the normTest function chosen by the user.}
}

\author{
Tomasz Zadlo}

\examples{
library(lme4)
mod <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
normCholTest(mod, shapiro.test)
} 


