\name{srswrRe}
\alias{srswrRe}
\title{Bootstrap sample of predicted random effects}
\description{The function draws at random a simple random sample with replacement from predicted random effects, where the sample size is equal the number of random effects in the whole population.
}
\usage{srswrRe(listRanef, reg)
}

\arguments{
  \item{listRanef}{{\emph{ranef(model)}} object where {\emph{model}} is an { \emph{lmer}} object.}
  \item{reg}{the population matrix of auxiliary variables named in \emph{fixed.part} and \emph{random.part}.}
}

\value{
\item{tablsrswrRe}{a vector of a simple random sample with replacement from predicted random effects, where the sample size is equal the number of random effects in the whole population.}
\item{lsrswrRe}{ a list of length equal the number of grouping variables taken into account in the random part of the model. Each list consists of 4 sublists: $raneftotal - a vector of a simple random sample with replacement from all predicted random effects under the cosidered grouping variable, $ranefname - a name of the grouping variable, $k - the number of random effects under the considered grouping variable, $df - a data frame of predicted random effects under the considered grouping variable, $dfsamp - a data frame of a simple random sample with replacement from predicted random effects under the considered grouping variable.}  
}

\references{
1. Carpenter, J.R., Goldstein, H. and Rasbash, J. (2003), A novel bootstrap procedure for assessing the relationship between class size and achievement. Journal of the Royal Statistical Society: Series C (Applied Statistics), 52, 431-443. \cr 

2. Chambers, R. and Chandra, H. (2013) A Random Effect Block Bootstrap for Clustered Data, Journal of Computational and Graphical Statistics, 22(2), 452-470. \cr

3. Thai, H.-T., Mentre, F., Holford, N.H., Veyrat-Follet, C. and Comets, E. (2013), A comparison of bootstrap approaches for estimating uncertainty of parameters in linear mixed-effects models. Pharmaceutical Statistics, 12, 129-140. \cr
}

\author{
Alicja Wolny-Dominiak, Tomasz Zadlo}

\examples{
library(lme4)
data(invData)
# data from one period are considered: 
invData2018 <- invData[invData$year == 2018,] 
attach(invData2018)
N <- nrow(invData2018) # population size
n <- 100 # sample size

set.seed(12345)
sampled_elements <- sample(N,n)
reg <- invData2018[, -which(names(invData2018) == 'investments')]

detach(invData2018)

invData2018sample <- invData2018[sampled_elements,]
attach(invData2018sample)
model <- lmer(investments ~ newly_registered + (1|NUTS2) + (1|NUTS4type))
srswrRe(ranef(model),reg)$tablsrswrRe
srswrRe(ranef(model),reg)$lsrswrRe

detach(invData2018sample)
}


