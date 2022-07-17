\name{bootResMis}
\alias{bootResMis}

\title{Residual bootstrap estimators of prediction accuracy under the misspecified model
}

\description{The function computes values of residual bootstrap estimators of RMSE and QAPE prediction accuracy measures of two predictors under the model assumed for one of them.
}

\usage{
bootResMis(predictorLMM, predictorLMMmis, B, p, correction)
}

\arguments{
  \item{predictorLMM}{plugInLMM object, the first predictor used to define the bootstrap model.
}
  \item{predictorLMMmis}{plugInLMM object, the second predictor.
}
  \item{B}{number of iterations in the bootstrap procedure.
}
  \item{p}{orders of quantiles in the QAPE.
}
  \item{correction}{logical. If TRUE, both bootstrapped random effects and random components are tranformed to avoid the problem of underdispersion of residual bootstrap distributions (see Details).
}
}


\details{Residual bootstrap considered by Carpener, Goldstein and Rasbash (2003), Chambers and Chandra (2013) and Thai et al. (2013) is used. We use model specification used in \emph{predictorLMM}. To generate one bootstrap realization of the population vector of the variable of interest: (i) from the sample vector of predicted random components the simple random sample with replacement of population size is drawn at random, (ii) from the vector of predicted random effects the simple random sample with replacement of size equal the number of random effects in the whole population is drawn at random. If {\emph{correction}} is {\emph{TRUE}}, then predicted random effects are transformed as described in Carpener, Goldstein and Rasbash (2003) in Section 3.2 and predicted random components as presented in Chambers and Chandra (2013) in Section 2.2. We use the MSE estimator defined as the mean of squared bootstrap errors considered by Rao and Molina (2015) p. 141 given by equation (6.2.22). The QAPE is a quantile of absolute prediction error which means that at least p100\% of realizations of absolute prediction errors are smaller or equal to QAPE. It is estimated as a quantile of absolute bootstrap errors as proposed by Zadlo (2017) in Section 2. The prediction accuracy of two predictors \emph{predictorLMM} and \emph{predictorLMMmis} is estimated under the model specified in \emph{predictorLMM}.}
\value{
\item{estQAPElmm}{estimated value/s of QAPE of \emph{predictorLMM} - number of rows is equal the number of orders of quantiles to be considered (declared in \emph{p}), number of columns is equal the number of predicted characteristics (declared in \emph{thetaFun}).}
\item{estRMSElmm}{estimated value/s of RMSE of \emph{predictorLMM} (more than one value is computed if in
\emph{thetaFun} more than one population characteristic is defined).}
\item{estQAPElmmMis}{estimated value/s of QAPE of \emph{predictorLMMmis} - number of rows is equal the number of orders of quantiles to be considered (declared in \emph{p}), number of columns is equal the number of predicted characteristics (declared in \emph{thetaFun}).}
\item{estRMSElmmMis}{estimated value/s of RMSE of \emph{predictorLMMmis} (more than one value is computed if in \emph{thetaFun} more than one population characteristic is defined).}
\item{predictorLMMSim}{bootstrapped values of \emph{predictorLMM}.}
\item{predictorLMMmisSim}{bootstrapped values of \emph{predictorLMMmis}.}
\item{thetaSim}{bootstrapped values of the predicted population or subpopulation characteristic/s.}
\item{Ysim}{simulated values of the (possibly tranformed) variable of interest.}
\item{errorLMM}{differences between bootstrapped values of \emph{predictorLMM} and bootstrapped values of the predicted characteristic/s.}
\item{errorLMMmis}{differences between bootstrapped values of \emph{predictorLMMmis} and bootstrapped values of the predicted characteristic/s.}
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
random.part.mis <- '(1|NUTS4type)'

reg <- invData2018[, -which(names(invData2018) == 'investments')]
weights <- rep(1,N) # homoscedastic random components

# Characteristics to be predicted:
# values of the variable for last two population elements  
thetaFun <- function(x) {x[c(379:380)]}

predictorLMM <- plugInLMM(YS, fixed.part, random.part, reg, con, weights, backTrans, thetaFun)
predictorLMM$thetaP

predictorLMMmis<-plugInLMM(YS, fixed.part, random.part.mis, reg, con, weights, backTrans, thetaFun)
predictorLMMmis$thetaP

set.seed(123456)
### Estimation of prediction accuracy
est_accuracy <- bootResMis(predictorLMM, predictorLMMmis, 10, c(0.5,0.8), correction = TRUE)

# Estimation of prediction RMSE of predictorLMM 
est_accuracy$estRMSElmm

# Estimation of prediction RMSE of predictorLMMmis
est_accuracy$estRMSElmmMis

# Estimation of prediction QAPE of predictorLMM 
est_accuracy$estQAPElmm

# Estimation of prediction QAPE of predictorLMMmis
est_accuracy$estQAPElmmMis

detach(invData2018)
}




