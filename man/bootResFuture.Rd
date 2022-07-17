\name{bootResFuture}
\alias{bootResFuture}

\title{
Residual bootstrap estimators of prediction accuracy - parallel computing
}

\description{The function computes values of residual bootstrap estimators of RMSE and QAPE prediction accuracy measures using parallel computing.
}

\usage{
bootResFuture(predictor, B, p, correction)
}

\arguments{
  \item{predictor}{one of objects: EBLUP, ebpLMMne or plugInLMM.
}
  \item{B}{number of iterations in the bootstrap procedure.
}
  \item{p}{orders of quantiles in the QAPE.
}
  \item{correction}{logical. If TRUE, both bootstrapped random effects and random components are tranformed to avoid the problem of underdispersion of residual bootstrap distributions (see Details).
}
}

\details{Residual bootstrap considered by Carpener, Goldstein and Rasbash (2003), Chambers and Chandra (2013) and Thai et al. (2013) is used. To generate one bootstrap realization of the population vector of the variable of interest: (i) from the sample vector of predicted random components the simple random sample with replacement of population size is drawn at random, (ii) from the vector of predicted random effects the simple random sample with replacement of size equal the number of random effects in the whole population is drawn at random. If {\emph{correction}} is {\emph{TRUE}}, then predicted random effects are transformed as described in Carpener, Goldstein and Rasbash (2003) in Section 3.2 and predicted random components as presented in Chambers and Chandra (2013) in Section 2.2. We use the MSE estimator defined as the mean of squared bootstrap errors considered by Rao and Molina (2015) p. 141 given by equation (6.2.22). The QAPE is a quantile of absolute prediction error which means that at least p100\% of realizations of absolute prediction errors are smaller or equal to QAPE. It is estimated as a quantile of absolute bootstrap errors as proposed by Zadlo (2017) in Section 2. The parallel processing is performed via the \emph{future.apply} package.}

\value{
\item{estQAPE}{estimated value/s of QAPE - number of rows is equal the number of orders of quantiles to be considered (declared in \emph{p}), number of columns is equal the number of predicted characteristics (declared in  in \emph{thetaFun}).}
\item{estRMSE}{estimated value/s of RMSE (more than one value is computed if in \emph{thetaFun} more than one population characteristic is defined).}
\item{predictorSim}{bootstrapped values of the predictor/s.}
\item{thetaSim}{bootstrapped values of the predicted population or subpopulation characteristic/s.}
\item{Ysim}{simulated values of the (possibly tranformed) variable of interest.}
\item{error}{differences between bootstrapped values of the predictor/s 
and bootstrapped values of the predicted characteristic/s.}
}

\references{
1. Carpenter, J.R., Goldstein, H. and Rasbash, J. (2003), A novel bootstrap procedure for assessing the relationship between class size and achievement. Journal of the Royal Statistical Society: Series C (Applied Statistics), 52, 431-443. \cr 

2. Chambers, R. and Chandra, H. (2013) A Random Effect Block Bootstrap for Clustered Data, Journal of Computational and Graphical Statistics, 22(2), 452-470. \cr

3. Thai, H.-T., Mentre, F., Holford, N.H., Veyrat-Follet, C. and Comets, E. (2013), A comparison of bootstrap approaches for estimating uncertainty of parameters in linear mixed-effects models. Pharmaceutical Statistics, 12, 129-140. \cr

}

\author{
Alicja Wolny-Dominiak, Tomasz Zadlo}

\examples{
data(invData) 
# data from one period are considered: 
invData2018 <- invData[invData$year == 2018,] 
attach(invData2018)

N <- nrow(invData2018) # population size

con <- rep(1,N) 
con[c(379:380)] <- 0 # last two population elements are not observed 

YS <- log(investments[con == 1]) # log-transformed values
backTrans <- function(x) exp(x) # back-transformation of the variable of interest
fixed.part <- 'log(newly_registered)'
random.part <- '(1|NUTS2)' 

reg <- invData2018[, -which(names(invData2018) == 'investments')]
weights <- rep(1,N) # homoscedastic random components

# Characteristics to be predicted:
# values of the variable for last two population elements  
thetaFun <- function(x) {x[c(379:380)]}
set.seed(123456)

predictor <- plugInLMM(YS, fixed.part, random.part, reg, con, weights, backTrans, thetaFun)
predictor$thetaP

### Estimation of prediction accuracy
est_accuracy <- bootResFuture(predictor, 10, c(0.5,0.8), correction = TRUE)

# Estimation of prediction RMSE
est_accuracy$estRMSE

# Estimation of prediction QAPE
est_accuracy$estQAPE

#         [,1]      [,2]
#50%  612.6089  67.45543
#80% 1886.9269 120.16246
####### Interpretations in case of prediction of investments 
####### for population element no. 379:
### It is estimated that at least 50% of absolute prediction errors are
# smaller or equal 612.6089 milion Polish zloty 
# and at least 50% of absolute prediction errors are
# greater or equal 612.6089 milion Polish zloty. 
### It is estimated that at least 80% of absolute prediction errors are
# smaller or equal 1886.9269 milion Polish zloty
# and at least 20% of absolute prediction errors are
# greater or equal 1886.9269 milion Polish zloty. 

detach(invData2018)
}


