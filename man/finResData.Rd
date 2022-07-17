\name{finResData}
\alias{finResData}
\docType{data}
\title{
Population data - financial results of enterprises in Poland at NUTS 2 level
}
\description{A data frame with 256 observations on the following 6 variables on 
Polish non-financial enterprises covering economic entities keeping accounting 
ledgers and which are obliged to prepare quarterly reports on income, costs 
and the financial result, with 50 and more persons.}

\arguments{
\item{year}{year.} 
\item{sector}{a sector (public or private) for which data are agregated.} 
\item{NUTS2}{NUTS 2 code (voivodships).}
\item{NUTS1}{NUTS 1 code (macroregions).}
\item{curAss}{current assets from period \emph{year-1} (real values).}
\item{netFinRes}{net financial result from period: \emph{year} (artificial variable).} 
}

\source{Statitics Poland, https://bdl.stat.gov.pl
}

\examples{

library(dplyr)
library(lme4)
library(qape) 
library(future.apply) 

data(finResData)
subpopulation <- interaction(finResData$sector, finResData$NUTS2)
finResData$subpopulation <- subpopulation

Ypop <- finResData$netFinRes

# netFinRes from years before 2021 are available in the sample: 
con <- ifelse(finResData$year != 2021, 1, 0)
YS <- Ypop[con == 1]

fixed.part <- 'curAss'
fixed.part.mis <- '1'
random.part <- '(1|subpopulation)'

reg <- select(finResData, -netFinRes)
weights <- rep(1, nrow(reg))
p <- c(0.5, 0.9)

# Prediction for first two subpopulations in 2021:
thetaFun <- function(x) {x[c(1,2)]}

#The predictor
plugin1 <- plugInLMM(YS, fixed.part, random.part, reg,
                     con, weights, backTrans = NULL, thetaFun)

}