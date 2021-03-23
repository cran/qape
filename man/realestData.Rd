\name{realestData}
\alias{realestData}
\docType{data}
\title{
Population data - real estate in Poland at NUTS 4 level
}
\description{A data frame with 1504 observations on the following 7 variables (NUTS 4 units with masked values of the variables due to Statistical confidentiality has been removed).}

\arguments{
\item{year}{year} 
\item{NUTS4}{NUTS 4 code (powiats)}
\item{NUTS2}{NUTS 2 code (voivodships)}
\item{NUTS4type}{type of NUTS 4 (1 - land counties, 2 - city counties/cities with powiat status)}
\item{premises}{number of residential premises sold in market transactions (in thousands)}
\item{area}{usable floor area of residential premises sold in market transactions (in millions of square meters)}
\item{price}{sum of prices of residential premises sold (in billions of Polish zloty)}
}

\source{Statitics Poland, https://bdl.stat.gov.pl/eng
}
\examples{
data(realestData)
hist(realestData$price[realestData$year==2018])
boxplot(realestData$price~realestData$year)
boxplot(realestData$price[realestData$year==2018]~realestData$NUTS2[realestData$year==2018])
boxplot(realestData$price[realestData$year==2018]~realestData$NUTS4type[realestData$year==2018])


library(lme4)
attach(realestData)

N <- nrow(realestData[(year == 2015),]) # population size in the first period
n <- 75 # sample size in the first period
set.seed(123456)
sampled_elements_in_2015 <- sample(N,n) 
con2015 <- rep(0,N)
con2015[sampled_elements_in_2015] <- 1

sampled_elements_in_2016 <- sample(N,n)
con2016 <- rep(0,N)
con2016[sampled_elements_in_2016] <- 1

sampled_elements_in_2017 <- sample(N,n)
con2017 <- rep(0,N)
con2017[sampled_elements_in_2017] <- 1

sampled_elements_in_2018 <- sample(N,n)
con2018 <- rep(0,N)
con2018[sampled_elements_in_2018] <- 1

con=as.logical(con2015, con2016, con2017, con2018)

model1 <- lmer(price ~ premises + area + (1|NUTS2)+(0+premises|NUTS2) + 
(1|NUTS4type)+(0+area|NUTS4type), subset=con)
AIC(model1)
model2 <- lmer(price ~ premises + area + (0+premises|NUTS2) + (0+area|NUTS4type), subset = con)
AIC(model2)

}