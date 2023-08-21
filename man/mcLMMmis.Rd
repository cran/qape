\name{mcLMMmis}
\alias{mcLMMmis}
\title{Monte Carlo simuation study of accuracy of predictors under the misspecified model
}
\description{The function computes in the Monte Carlo simulation study values of accuracy measures of three predictors under the model assumed for one of them with possible modifications of covariance matrices of random effects and random components.
}
\usage{
mcLMMmis(Ypop, predictorLMMmis, predictorLMM, predictorLMM2, K, p, ratioR, ratioG)
}
\arguments{
  \item{Ypop}{population values of the variable of interest (already transformed if necessary) which are used as the dependent variable in the population model.
}
  \item{predictorLMMmis}{plugInLMM object, the predictor used to define the model assumed in the simulation study.
}
  \item{predictorLMM}{plugInLMM object, the first predictor, the accuracy of which is assessed in the simulation study.
}
  \item{predictorLMM2}{plugInLMM object, the second predictor, the accuracy of which is assessed in the simulation study.
}
  \item{K}{the number of Monte Carlo iterations.
}
  \item{p}{orders of quantiles in the QAPE.
}
   \item{ratioR}{the value by which the diagonal elements of the covariance matrix of random components of the model based on the whole population data and formulation used in \emph{predictorLMMmis} are divided. Then, the corrected covariance matrix is used to generate bootstrap realizations of random components.
}
  \item{ratioG}{the value by which the diagonal elements of the covariance matrix of random effects of the model based on the whole population data and formulation used in \emph{predictorLMMmis} are divided. Then, the corrected covariance matrix, assuming that it is positive definite, is used to generate bootstrap realizations of random effects. If it is not positive definite, the alert is printed and the dependent variable is generated based on the model without random effects.
}
}

\details{In the model-based simulation study population values of the dependent variable are generated based on the (possibly transformed) Linear Mixed Model used in \emph{predictorLMMmis} with possibly modified covariance matrices of random effects and random components by the usage of \emph{ratioR} and \emph{ratioG} arguments. In the simulation study accuracy of predictors \emph{predictorLMM} and \emph{predictorLMM2} is assessed. Although, all the predictors are \emph{plugInLMM} objects, it should be noted that under the non-transformed Linear Mixed Model and in the case of the prediction of the linear combination of the dependent variable (e.g. the mean, the total, and one realization of the variable), the predictors are Empirical Best Linear Unbiased Predictors. What is more, if \emph{predictorLMMmis} is defined as \emph{predictorLMM}, the Monte Carlo simulation study of accuracy of \emph{predictorLMM} under correctly specified model and of \emph{predictorLMM2} under misspecified model is conducted.      
}


\value{
\item{errorLMM}{Monte Carlo prediction errors of \emph{predictorLMM} - number of rows is equal to the number of predicted characteristics (declared in \emph{thetaFun}), number of columns is equal to \emph{K}.} 

\item{errorLMM2}{Monte Carlo prediction errors of \emph{predictorLMM2} - number of rows is equal to the number of predicted characteristics (declared in \emph{thetaFun}), number of columns is equal to \emph{K}.} 

\item{QAPElmm}{value/s of the QAPE of \emph{predictorLMM} assessed in the Monte Carlo study - number of rows is equal the number of orders of quantiles to be considered (declared in \emph{p}), number of columns is equal the number of predicted characteristics (declared in \emph{thetaFun}).}

\item{RMSElmm}{value/s of the RMSE of \emph{predictorLMM} assessed in the Monte Carlo study (more than one value is computed if in \emph{thetaFun} more than one population characteristic is defined).}

\item{rRMSElmm}{value/s of the rRMSE (in percentages) of \emph{predictorLMM} assessed in the Monte Carlo study (more than one value is computed if in \emph{thetaFun} more than one population characteristic is defined).}

\item{rBlmm}{value/s of the relative bias (in percentages) of \emph{predictorLMM} assessed in the Monte Carlo study (more than one value is computed if in \emph{thetaFun} more than one population characteristic is defined).}

\item{QAPElmm2}{value/s of the QAPE of \emph{predictorLMM2} assessed in the Monte Carlo study - number of rows is equal the number of orders of quantiles to be considered (declared in \emph{p}), number of columns is equal the number of predicted characteristics (declared in \emph{thetaFun}).}

\item{RMSElmm2}{value/s of the RMSE of \emph{predictorLMM2} assessed in the Monte Carlo study (more than one value is computed if in \emph{thetaFun} more than one population characteristic is defined).}

\item{rRMSElmm2}{value/s of the rRMSE (in percentages) of \emph{predictorLMM2} assessed in the Monte Carlo study (more than one value is computed if in \emph{thetaFun} more than one population characteristic is defined).}

\item{rBlmm2}{value/s of the relative bias (in percentages) of \emph{predictorLMM2} assessed in the Monte Carlo study (more than one value is computed if in \emph{thetaFun} more than one population characteristic is defined).}

\item{positiveDefiniteEstG}{logical indicating if the estimated covariance matrix of random effects, used to generate Monte Carlo realizations of the dependent variable, is positive definite.}
}


\references{

1. Chatterjee, S., Lahiri, P. Li, H. (2008) Parametric bootstrap approximation to the distribution of EBLUP and related prediction intervals in linear mixed models, Annals of Statistics, Vol. 36 (3), 1221-1245.\cr

2. Rao, J.N.K. and Molina, I. (2015) Small Area Estimation. Second edition, John Wiley & Sons, New Jersey.\cr

3. Zadlo T. (2017), On asymmetry of prediction errors in small area estimation, Statistics in Transition, 18 (3), 413-432. \cr

}

\author{
Tomasz Zadlo}

\examples{

data(invData) 
# data from one period are considered: 
invData2018 <- invData[invData$year == 2018,] 

invData2018$investments <- invData2018$investments/1000
attach(invData2018)
N <- nrow(invData2018) # population size

con <- rep(0,N) 
set.seed(123456)
con[sample(N,50)] <- 1 # sample size equals 50

YS <- log((investments[con == 1])) # log-transformed values
backTrans <- function(x) {exp(x)} # back-transformation of the variable of interest
fixed.part <- 'log(newly_registered)'
random.part <- '(1|NUTS2)'
random.part2 <- '(1|NUTS4type)'

reg <- invData2018[, -which(names(invData2018) == 'investments')]
weights <- rep(1,N) # homoscedastic random components
weights.mis <- sqrt(newly_registered)

# Characteristics to be predicted:
# the population mean and the population total  
thetaFun <- function(x) {c(mean(x), median(x))}

predictorLMMmis <- plugInLMM(YS, fixed.part, random.part, reg,con,weights.mis,backTrans,thetaFun)
predictorLMMmis$thetaP

predictorLMM <- plugInLMM(YS, fixed.part, random.part, reg, con, weights, backTrans, thetaFun)
predictorLMM$thetaP

predictorLMM2 <- plugInLMM(YS, fixed.part, random.part2, reg, con, weights, backTrans, thetaFun)
predictorLMM2$thetaP

Ypop <- log(invData2018$investments)

# Monte Carlo simulation study under the misspecified model defined in predictorLMMmis 
# with modified covariance matrices R and G
set.seed(123456)
mcLMMmis(Ypop, predictorLMMmis, predictorLMM, predictorLMM2, 5, c(0.75,0.9), 2, 0.1)

# Monte Carlo simulation study under the model defined in predictorLMM
# correctly specified for predictorLMM and misspecified for predictorLMM2
set.seed(123456)
mcLMMmis(Ypop, predictorLMM, predictorLMM, predictorLMM2, 5, c(0.75,0.9), 1, 1)

detach(invData2018)
}




