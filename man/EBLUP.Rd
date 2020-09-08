\name{EBLUP}
\alias{EBLUP}
\title{
Empirical Best Linear Unbiased Predictor
}
\description{
The function computes the value of the EBLUP of the linear combination of the variable of interest under linear mixed model estimated using REML.
}
\usage{
EBLUP(YS, fixed.part, random.part, reg, con, gamma, weights)
}

\arguments{
  \item{YS}{values of the variable of interest observed in the sample.
}
  \item{fixed.part}{fixed-effects terms declared as in \emph{lmer} object. 
}
  \item{random.part}{random-effects terms declared as in \emph{lmer} object. 
}
  \item{reg}{the population matrix of auxiliary variables named in \emph{fixed.part} and \emph{random.part}.
}
  \item{con}{the population 0-1 vector with 1s for elements in the sample and 0s for elements which are not in the sample.
}
  \item{gamma}{the population vector which transpose multiplied by the population vector of the variable of interest gives the predicted characteristic. For example, if \emph{gamma} is the population vector on 1s, the sum of the values of the variable of interest in the whole dataset is predicted.  
}
  \item{weights}{the population vector of weights, defined as in \emph{lmer} object, allowing to include heteroscedasticity of random components in the mixed linear model.  
}
}

\details{The function computes the value of the EBLUP of the linear combination of the variable of interest based on the formula (21) in Zadlo (2017) (see Remark 5.1 in the paper for further explanations). Predicted values for unsampled population elements in subsets for which random effects are not observed in the sample are computed based only on fixed effects.
}
\value{
The function returns a list with the following objects:
	\item{fixed.part}{the fixed part of the formula of model.}
        \item{random.part}{the random part of the formula of model.} 
        \item{thetaP}{the value of the predictor.}
        \item{beta}{the estimated vector of fixed effects.}
        \item{Xbeta}{the product of two matrices: the population model matrix of auxiliary variables X and the estimated vector of fixed effects.}
        \item{sigma2R}{the estimated variance parameter of the distribution of random components.}
        \item{R}{the estimated covariance matrix of random components for sampled elements.}
        \item{G}{the estimated covariance matrix of random effects.}
        \item{model}{the formula of the model (as in \emph{lmer} object).}
        \item{mEst}{\emph{lmer} object with the estimated model.}
        \item{YS}{the sample vector of the variable of interest.}
        \item{reg}{the population matrix of auxiliary variables named in \emph{fixed.part} and \emph{random.part}.}
        \item{con}{the population 0-1 vector with 1s for elements in the sample and 0s for elements which are not in the sample.}
        \item{regS}{the sample matrix of auxiliary variables named in \emph{fixed.part} and \emph{random.part}.}
        \item{regR}{the matrix of auxiliary variables named in \emph{fixed.part} and \emph{random.part} for population elements which are not observed in the sample.}
        \item{gamma}{the population vector which transpose multiplied by the population vector of the variable of interest gives the predicted characteristic.}
        \item{gammaS}{the subvector of \emph{gamma} for sampled elements.}
        \item{gammaR}{the subvector of \emph{gamma} for population elements which are not observed in the sample.}
        \item{weights}{the population vector of weigts, defined as in \emph{lmer} object, allowing to include the heteroscedasticity of random components in the mixed linear model.}
        \item{Z}{the population model matrix of auxiliary variables associated with random effects.}
        \item{X}{the population model matrix of auxiliary variables associated with fixed effects.}
        \item{ZS}{the submatrix of Z matrix where the number of rows equals the number of sampled elements and the number of columns equals the number of estimated random effects.}
        \item{XR}{the submatrix of X matrix (with the same number of columns) for population elements which are not observed in the sample.}
        \item{ZR}{the submatrix of Z matrix where the number of rows equals the number of population elements which are not observed in the sample and the number of columns equals the number of estimated random effects.}
        \item{eS}{the sample vector of estimated random components.}
        \item{vS}{the estimated vector of random effects.}
}
\references{
1. Henderson, C.R. (1950) Estimation of Genetic Parameters (Abstract). Annals of Mathematical Statistics 21, 309-310. \cr
2. Royall, R.M. (1976) The Linear Least Squares Prediction Approach to Two-Stage Sampling. Journal of the American Statistical Association 71, 657-473. \cr
3. Zadlo, T. (2017) On prediction of population and subpopulation characteristics for future periods, Communications in Statistics - Simulation and Computation 461(10), 8086-8104. \cr
}
\author{
Alicja Wolny-Dominiak, Tomasz Zadlo}

\examples{

library(lme4)
library(Matrix)


### Prediction of the subpopulation mean based on the cross-sectional data

data(invData) 
#data from one period are considered: 
invData2018<-invData[invData$year == 2018,] 
attach(invData2018)

N=nrow(invData2018) #population size
n=100 # sample size
#subpopulation of interest: NUTS4type==2
Nd=sum(NUTS4type==2) #subpopulation size

set.seed(123456)
sampled_elements=sample(N,n)
con=rep(0,N)
con[sampled_elements]=1 #elements in the sample
YS=investments[sampled_elements]
fixed.part <- 'newly_registered'
random.part <- '(newly_registered| NUTS2)' 
reg=invData2018[, - which(names(invData2018) == 'investments')]

gamma=rep(0,N)
gamma[NUTS4type==2]=1/Nd

weights=rep(1,N) #homoscedastic random components

#Predicted value of the mean in the following subpopulation: NUTS4type==2
EBLUP(YS, fixed.part, random.part, reg, con, gamma, weights)$thetaP

#All results
EBLUP(YS, fixed.part, random.part, reg, con, gamma, weights)

detach(invData2018)

##########################################################

### Prediction of the subpopulation total based on longitudinal data

data(invData)
attach(invData)

N=nrow(invData[(year==2013),]) #population size in the first period
n=38 # sample size in the first period
#subpopulation and time period of interest: NUTS2=='02' & year==2018
Ndt=sum(NUTS2=='02' & year==2018) #subpopulation size in the period of interest

set.seed(123456)
sampled_elements_in_2013=sample(N,n)
con2013=rep(0,N)
con2013[sampled_elements_in_2013]=1 #elements in the sample in 2013

#balanced panel sample - the same elements in all 6 periods:
con=rep(con2013,6)

YS=investments[con==1]
fixed.part <- 'newly_registered'
random.part <- '(newly_registered  | NUTS4)' 
reg=invData[, - which(names(invData) == 'investments')]

gamma=rep(0,nrow(invData))
gamma[NUTS2=='02' & year==2018]=1

weights=rep(1,nrow(invData)) #homoscedastic random components

# Predicted value of the total 
# in the following subpopulation: NUTS4type==2 
# in the following time period: year==2018
EBLUP(YS, fixed.part, random.part, reg, con, gamma, weights)$thetaP

#All results
EBLUP(YS, fixed.part, random.part, reg, con, gamma, weights)

detach(invData)

}