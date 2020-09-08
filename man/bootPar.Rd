\name{bootPar}
\alias{bootPar}
\title{
Parametric bootstrap estimators of prediction accuracy
}
\description{The function computes values of parametric bootstrap estimators of RMSE and QAPE prediction accuracy measures.
}
\usage{
bootPar(predictor, B, p)
}

\arguments{
  \item{predictor}{one of objects: EBLUP, ebpLMMne or plugInLMM.
}
  \item{B}{number of iterations in the bootstrap procedure.
}
  \item{p}{orders of quantiles in the QAPE.
}
}

\details{We use bootstrap model presented by Chatterjee, Lahiri and  Li (2008) p. 1229 but assumed for all population elements. Vectors of random effects and random components are generated from the multivariate normal distribution where REML estimates of model parameters are used. Random effects are generated for all population elements even for subsets with zero sample sizes (for which random effects are not estimated). We use the MSE estimator defined as the mean of squared bootstrap errors considered by Rao and MOlina (2015) p. 141 and given by equation (6.2.22). The QAPE is a quantile of absolute prediction error which means that at least p100\% of realizations of absolute prediction errors are smaller or equal to QAPE. It is estimated as a quantile of absolute bootstrap errors as proposed by Zadlo (2017) in Section 2.  
}

\value{
\item{estQAPE}{estimated value/s of QAPE - number of rows is equal the number of orders of quantiles to be considered (declared in \emph{p}), number of columns is equal the number of predicted characteristics (declared in  in \emph{thetaFun}).}
\item{estRMSE}{estimated value/s of RMSE (more than one value is computed if in \emph{thetaFun} more than one population characteristic is defined).}
\item{predictorSim}{bootstrapped values of the predictor/s.}
\item{thetaSim}{bootstrapped values of the predicted population or subpopulation characteristic/s.}
\item{YSim}{simulated values of the (possibly tranformed) variable of interest.}
\item{error}{differences between bootstrapped values of the predictor/s and bootstrapped values of the predicted characteristic/s.}
}

\references{

1. Butar, B. F., Lahiri, P. (2003) On measures of uncertainty of empirical Bayes small-area estimators, Journal of Statistical Planning and Inference, Vol. 112, pp. 63-76.\cr

2. Chatterjee, S., Lahiri, P. Li, H. (2008) Parametric bootstrap approximation to the distribution of EBLUP and related prediction intervals in linear mixed models, Annals of Statistics, Vol. 36 (3), pp. 1221-1245.\cr

3. Rao, J.N.K. and Molina, I. (2015) Small Area Estimation. Second edition, John Wiley & Sons, New Jersey.\cr

4. Zadlo T. (2017), On asymmetry of prediction errors in small area estimation, Statistics in Transition, 18 (3), 413-432 \cr

}

\author{
Alicja Wolny-Dominiak, Tomasz Zadlo}

\examples{

library(lme4)
library(Matrix)
library(mvtnorm)


data(invData) 
#data from one period are considered: 
invData2018<-invData[invData$year == 2018,] 
attach(invData2018)

N=nrow(invData2018) #population size

con=rep(1,N) 
con[c(379,380)]<-0 # last two population elements are not observed 

YS=(investments[con==1]) # log-transformed values
backTrans <- function(x) x # back-transformation of the variable of interest
fixed.part <- 'log(newly_registered)'
random.part <- '((1|NUTS2)+((newly_registered-1)|NUTS2))' 

reg=invData2018[, - which(names(invData2018) == 'investments')]
weights=rep(1,N) #homoscedastic random components

# Characteristics to be predicted:
# values of the variable for last two population elements  
thetaFun <- function(x) {x[c(379,380)]}
set.seed(123456)

# Predicted values of quartiles 
# in the following subpopulation: NUTS4type==2 
# in the following time period: year==2018

predictor=plugInLMM(YS, fixed.part, random.part, reg, con, weights, backTrans, thetaFun)
predictor$thetaP

# Estimation of prediction accuracy
est_accuracy=bootPar(predictor, 10, c(0.75,0.9))

# Estimation of prediction RMSE
est_accuracy$estRMSE

# Estimation of prediction QAPE
est_accuracy$estQAPE

}

