\name{invData}
\alias{invData}
\docType{data}
\title{
Population data - investments in Poland at NUTS 4 level
}
\description{A data frame with 2280 observations on 6 variables presented below.}

\arguments{
\item{year}{year.} 
\item{NUTS4}{NUTS 4 code (powiats).}
\item{NUTS2}{NUTS 2 code (voivodships).}
\item{NUTS4type}{type of NUTS 4 (1 - land counties, 2 - city counties/cities with powiat status).}
\item{investments}{investment outlays in millions PLN, in current prices; data concern Polish economic entities, including independent health care facilities and cultural institutions with legal personalities in which the number of employed persons exceeds 9 (source of data: Annual survey of the economic activity of enterprises conducted by Statistics Poland).}
\item{newly_registered}{newly registered entities of the national economy recorded in the REGON register (in thousands).}
}

\source{Statistics Poland, https://bdl.stat.gov.pl/eng
}
\examples{
data(invData)
hist(invData$newly_registered[invData$year==2018])
boxplot(invData$investments~invData$year)
boxplot(invData$investments[invData$year==2018]~invData$NUTS2[invData$year==2018])
boxplot(invData$investments[invData$year==2018]~invData$NUTS4type[invData$year==2018])
}