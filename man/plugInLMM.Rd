\name{plugInLMM}
\alias{plugInLMM}
\title{
PLUG-IN predictor based on the linear mixed model
}
\description{
The function computes the value of the plug-in predictor under the linear mixed model estimated using REML assumed for possibly transformed variable of interest. 
}
\usage{
plugInLMM(YS, fixed.part, random.part, reg, con, weights, backTrans, thetaFun)
}

\arguments{
  \item{YS}{values of the variable of interest (already transformed if necessary) observed in the sample and used in the model as the dependent variable.}
  
  \item{fixed.part}{fixed-effects terms declared as in \emph{lmer} object.
}
  \item{random.part}{random-effects terms declared as in \emph{lmer} object.
}
  \item{reg}{the population matrix of auxiliary variables named in \emph{fixed.part} and \emph{random.part}.
}
  \item{con}{the population 0-1 vector with 1s for elements in the sample and 0s for elements which are not in the sample.
}
  \item{weights}{the population vector of weights, defined as in \emph{lmer} object, allowing to include the heteroscedasticity of random components in the mixed linear model.
}
  \item{backTrans}{back-transformation function of the variable of interest (e.g. if YS is log-tranformed, then backTrans <- function(x) exp(x)). 
}
  \item{thetaFun}{the predictor function (e.g. mean or sd)
}
}

\details{The function computes the value of the plug-in estimator in two steps as presented by Chwila and Zadlo (2019) p. 20. Firstly, we build the population vector consisting of real values of the variable of interest for sampled elements and (possibly back-transformed) fitted values of the variable of interest based on the estimated model. Secondly, the value/s of {\emph{thetaFun}} based on the population vector built in the first step is/are computed. Predicted values for unsampled population elements in subsets for which random effects are not observed in the sample are computed based only on fixed effects.
}

\value{
The function returns a list with the following objects:
\item{thetaP}{the value/s of the predictor (more than one value is computed if in \emph{thetaFun} more than one population characteristic is defined).}
\item{fixed.part}{the fixed part of the formula of model.}
\item{random.part}{the random part of the formula of model.} 
\item{thetaFun}{the function of the population values of the variable of interest (on the original scale) which defines at least one population or subpopulation characteristic to be predicted.}
\item{backTrans}{back-transformation function of the variable of interest (e.g. if YS used in the model is log-tranformed, then backTrans <- function(x) exp(x)).}
\item{YP}{predicted values of the variable of interest for unsampled elements (without back-tranformation).}
\item{YbackTrans}{population vector of the values of the variable of interest on the orignal scale for sampled elements and back-transformed predicted values of the variable of interest for unsampled elements.}
\item{YPbackTrans}{back-transformed predicted values of the variable of interest for unsampled elements.}
\item{beta}{the estimated vector of fixed effects.}
\item{Xbeta}{the product of two matrices: the population model matrix of auxiliary variables X and the estimated vector of fixed effects.}
\item{sigma2R}{the estimated variance parameter of the distribution of random components.}
\item{R}{the estimated covariance matrix of random components for sampled elements.}
\item{G}{the estimated covariance matrix of random effects.}
\item{model}{the formula of the model (as in \emph{lmer} object).}
\item{mEst}{\emph{lmer} object with the estimated model.}
\item{YS}{values of the variable of interest (already transformed if necessary) observed in the sample and used in the model as the dependent variable.}
\item{reg}{the population matrix of auxiliary variables named in \emph{fixed.part} and \emph{random.part}.}
\item{con}{the population 0-1 vector with 1s for elements in the sample and 0s for elements which are not in the sample.}
\item{regS}{the sample matrix of auxiliary variables named in \emph{fixed.part} and \emph{random.part}.}
\item{regR}{the matrix of auxiliary variables named in \emph{fixed.part} and \emph{random.part} for unsampled population elements.}
\item{weights}{the population vector of weigts, defined as in \emph{lmer} object, allowing to include the heteroscedasticity of random components in the mixed linear model.}
\item{Z}{the population model matrix of auxiliary variables associated with random effects.}     \item{ZBlockNames}{labels of blocks of random effects in Z matrix.}
\item{ZS}{the submatrix of Z matrix where the number of rows equals the number of sampled elements and the number of columns equals the number of estimated random effects.}
\item{XR}{the submatrix of X matrix (with the same number of columns) for unsampled population elements.}
\item{ZR}{the submatrix of Z matrix where the number of rows equals the number of unsampled population elements and the number of columns equals the number of estimated random effects.}
\item{eS}{the sample vector of estimated random components.}
\item{vS}{the estimated vector of random effects.}
}

\references{
Chwila, A., Zadlo, T. (2019) On properties of empirical best predictors. Communications in Statistics - Simulation and Computation, available online: https://doi.org/10.1080/03610918.2019.1649422

}
\author{
Alicja Wolny-Dominiak, Tomasz Zadlo}

\examples{

library(lme4)
library(Matrix)


### Prediction of the subpopulation median 
### and the subpopulation standard deviation 
### based on the cross-sectional data

data(invData) 
# data from one period are considered: 
invData2018 <- invData[invData$year == 2018,] 
attach(invData2018)

N <- nrow(invData2018) # population size
n <- 100 # sample size

set.seed(123456)
sampled_elements <- sample(N,n)
con <- rep(0,N)
con[sampled_elements] <- 1 # elements in the sample
YS <- log(investments[sampled_elements]) # log-transformed values
backTrans <- function(x) exp(x) # back-transformation of the variable of interest
fixed.part <- 'log(newly_registered)'
random.part <- '(log(newly_registered)  | NUTS2)' 
reg <- invData2018[, -which(names(invData2018) == 'investments')]
weights <- rep(1,N) # homoscedastic random components

# Characteristics to be predicted - the median and the standard deviation
# in following subpopulation: NUTS4type == 2
thetaFun <- function(x) {c(median(x[NUTS4type == 2]),sd(x[NUTS4type == 2]))}

# Predicted values of the median and the standard deviation
# in the following subpopulation: NUTS4type == 2

plugInLMM(YS, fixed.part, random.part, reg, con, weights, backTrans, thetaFun)$thetaP

# All results
plugInLMM(YS, fixed.part, random.part, reg, con, weights, backTrans, thetaFun)

detach(invData2018)

##########################################################

### Prediction of the subpopulation quartiles based on longitudinal data

data(invData)
attach(invData)

N <- nrow(invData[(year == 2013),]) # population size in the first period
n <- 38 # sample size in the first period
# subpopulation and time period of interest: NUTS2 == '02' & year == 2018
Ndt=sum(NUTS2=='02' & year==2018) # subpopulation size in the period of interest

set.seed(123456)
sampled_elements_in_2013 <- sample(N,n)
con2013 <- rep(0,N)
con2013[sampled_elements_in_2013] <- 1 # elements in the sample in 2013

# balanced panel sample - the same elements in all 6 periods:
con <- rep(con2013,6)

YS <- log(investments[con == 1]) # log-transformed values
backTrans <- function(x) exp(x) # back-transformation of the variable of interest
fixed.part <- 'log(newly_registered)'
random.part <- '(0 + log(newly_registered) | NUTS4)' 
reg <- invData[, -which(names(invData) == 'investments')]
weights <- rep(1,nrow(invData)) # homoscedastic random components

# Characteristics to be predicted - quartiles in 2018 
# in the following subpopulation: NUTS4type == 2
thetaFun <- function(x) {quantile(x[NUTS2 == '02' & year == 2018],probs = c(0.25,0.5,0.75))}

# Predicted values of quartiles 
# in the following subpopulation: NUTS4type == 2 
# in the following time period: year == 2018
plugInLMM(YS, fixed.part, random.part, reg, con, weights, backTrans, thetaFun)$thetaP

# All results
plugInLMM(YS, fixed.part, random.part, reg, con, weights, backTrans, thetaFun)

detach(invData)

}


