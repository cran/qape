\name{bootParFutureCor}
\alias{bootParFutureCor}

\title{Parametric bootstrap estimators of prediction accuracy - parallel computing using corrected covariance matrices
}


\description{The function computes values of parametric bootstrap estimators of RMSE and QAPE prediction accuracy measures using parallel computing under the misspecified model. The model misspecification is obtained by the modification of the covariance matrices of random effects and random components estimated based on sample data. The correction is made by the division of the diagonal elements of random effects and random components estimated based on sample data by values defined by users and then, the corrected covariance matrices are used to generate bootstrap realizations of the dependent variables.   
}


\usage{
bootParFutureCor(predictor, B, p, ratioR, ratioG)
}

\arguments{
  \item{predictor}{one of objects: EBLUP, ebpLMMne or plugInLMM.
}
  \item{B}{number of iterations in the bootstrap procedure.
}

  \item{p}{orders of quantiles in the QAPE.
}


  \item{ratioR}{the value by which the diagonal elements of the covariance matrix of random components estimated based on sample data are divided. Then, the corrected covariance matrix is used to generate bootstrap realizations of random components.
}
  \item{ratioG}{the value by which the diagonal elements of the covariance matrix of random effects estimated based on sample data are divided. Then, the corrected covariance matrix, assuming that it is positive definite, is used to generate bootstrap realizations of random effects. If it is not positive definite, the alert is printed and the dependent variable is generated based on the model without random effects.
}
}

\details{We use bootstrap model presented by Chatterjee, Lahiri and  Li (2008) p. 1229 but assumed for all population elements. Vectors of random effects and random components are generated from the multivariate normal distribution, where REML estimates of model parameters are used. Random effects are generated for all population elements, even for subsets with zero sample sizes (for which random effects are not estimated). We use the MSE estimator defined as the mean of squared bootstrap errors considered by Rao and Molina (2015) p. 141 and given by equation (6.2.22). The QAPE is a quantile of absolute prediction error, which means that at least p100\% of realizations of absolute prediction errors are smaller or equal to QAPE. It is estimated as a quantile of absolute bootstrap errors, as proposed by Zadlo (2017) in Section 2.  The parallel processing is performed via the \emph{future.apply} package. The dependent variable is generated based on the modified (misspecified) model with corrected covariance matrices of random effects and random components. The correction is made by the division of the diagonal elements of the covariance matrix of random components estimated based on sample data by \emph{ratioR}, and by the division of the diagonal elements of the covariance matrix of random effects estimated based on sample data by \emph{ratioG}. If the estimated covariance matrix of random effect after the correction is not positive definite, the alert is printed and the bootstrap realizations of dependent variable are generated based on the model without random effects.    
}

\value{
\item{estQAPE}{estimated value/s of QAPE - number of rows is equal to the number of orders of quantiles to be considered (declared in \emph{p}), number of columns is equal to the number of predicted characteristics (declared in \emph{thetaFun}).}
\item{estRMSE}{estimated value/s of RMSE (more than one value is computed if in \emph{thetaFun} more than one population characteristic is defined).}
\item{predictorSim}{bootstrapped values of the predictor/s.}
\item{thetaSim}{bootstrapped values of the predicted population or subpopulation characteristic/s.}
\item{Ysim}{simulated values of the (possibly tranformed) variable of interest.}
\item{error}{differences between bootstrapped values of the predictor/s and bootstrapped values of the predicted characteristic/s.}
\item{positiveDefiniteEstG}{logical indicating if the estimated covariance matrix of random effects, used to generate bootstrap realizations of the dependent variable, is positive definite.}
}

\references{

1. Butar, B. F., Lahiri, P. (2003) On measures of uncertainty of empirical Bayes small-area estimators, Journal of Statistical Planning and Inference, Vol. 112, pp. 63-76.\cr

2. Chatterjee, S., Lahiri, P. Li, H. (2008) Parametric bootstrap approximation to the distribution of EBLUP and related prediction intervals in linear mixed models, Annals of Statistics, Vol. 36 (3), pp. 1221?1245.\cr

3. Rao, J.N.K. and Molina, I. (2015) Small Area Estimation. Second edition, John Wiley & Sons, New Jersey.\cr

4. Zadlo T. (2017), On asymmetry of prediction errors in small area estimation, Statistics in Transition, 18 (3), 413-432. \cr

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
con[c(379,380)] <- 0 # last two population elements are not observed 

YS <- log(investments[con == 1]) # log-transformed values
backTrans <- function(x) exp(x) # back-transformation of the variable of interest
fixed.part <- 'log(newly_registered)'
random.part <- '(1|NUTS2)'

reg <- invData2018[, -which(names(invData2018) == 'investments')]
weights <- rep(1,N) # homoscedastic random components

# Characteristics to be predicted:
# values of the variable for last two population elements  
thetaFun <- function(x) {x[c(379,380)]}
set.seed(123)

predictor <- plugInLMM(YS, fixed.part, random.part, reg, con, weights, backTrans, thetaFun)
predictor$thetaP

### Estimation of prediction accuracy under the misspecified model
est_accuracy <- bootParFutureCor(predictor, 10, c(0.75,0.9), 2, 0.01)

# Estimation of prediction RMSE under the misspecified model
est_accuracy$estRMSE

# Estimation of prediction QAPE under the misspecified model
est_accuracy$estQAPE

detach(invData2018)
}



