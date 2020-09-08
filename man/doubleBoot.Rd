\name{doubleBoot}
\alias{doubleBoot}

\title{Double bootstrap estimators of prediction accuracy}

\description{The function computes values of double bootstrap estimators of MSE and QAPE prediction accuracy measures.
}
\usage{
doubleBoot(predictor, B1, B2, p, q)
}

\arguments{
  \item{predictor}{one of objects: EBLUP, ebpLMMne or plugInLMM.
}
  \item{B1}{number of first-level bootstrap iterations.

}
  \item{B2}{
number of second-level bootstrap iterations.
}
  \item{p}{
orders of quantiles in the QAPE.
}
  \item{q}{
estimator bounds assumed for \emph{estMSE_db_1_EF} and \emph{estMSE_db_telesc_EF} (which are corrected versions of \emph{estMSE_db_1} and \emph{estMSE_db_telesc})
}
}

\details{Double-bootstrap method considered by Hall and Maiti (2006) and  Erciulescu and Fuller (2013) is used. Vectors of random effects and random components are generated from the multivariate normal distribution and REML estimates of model parameters are used. Random effects are generated for all population elements even for subsets with zero sample sizes (for which random effects are not estimated). Double-bootstrap MSE estimator presented in Hall and Maiti (2006) and  Erciulescu and Fuller (2013) are taken into account. The QAPE is a quantile of absolute prediction error which means that at least p100\% of realizations of absolute prediction errors are smaller or equal to QAPE.  
}


\value{

  \item{estMSE_param}{value/s of the parametric bootstrap MSE estimator. More than one value is computed if in \emph{thetaFun} more than one population characteristic is defined.
    }
  \item{estMSE_db_B2}{value/s of the double bootstrap MSE estimator computed as the difference of doubled value of \emph{estMSE_param} and the second-level MSE estimator based on B2 iterations. More than one value is computed if in \emph{thetaFun} more than one population characteristic is defined.
    }
  \item{estMSE_db_B2_WDZ}{value/s of the double bootstrap MSE estimator computed as the mean of squared first-level boostraped errors, each corrected by the mean of squared second-level boostraped errors based on B2 iterations (where correction is made only if their difference is non-negative). More than one value is computed if in \emph{thetaFun} more than one population characteristic is defined.
    }
  \item{estMSE_db_B2_HM}{value/s of the double bootstrap MSE estimator proposed by Hall and Maiti (2006) equation (2.17). More than one value is computed if in \emph{thetaFun} more than one population characteristic is defined.
    }
  \item{estMSE_db_1}{value/s of the double bootstrap MSE estimator computed as the difference of doubled value of \emph{estMSE_param} and the second-level MSE estimator based on B2=1 iteration. More than one value is computed if in \emph{thetaFun} more than one population characteristic is defined.
    }
  
  \item{estMSE_db_1_WDZ}{value/s of the double bootstrap MSE estimator computed as the mean of squared first-level boostraped errors, each corrected by the squared second-level boostraped error based on 1 iteration (where correction is made only if their difference is non-negative). More than one value is computed if in \emph{thetaFun} more than one population characteristic is defined.
    }
  \item{estMSE_db_1_EF}{value/s of the double bootstrap MSE estimator proposed by Erciulescu and Fuller (2014) given by equation (13) with correction (17), where the bound for the correction is declared as \emph{q}. More than one value is computed if in \emph{thetaFun} more than one population characteristic is defined.
    }
  \item{estMSE_db_telesc}{value/s of the telescoping double bootstrap MSE estimator proposed by Erciulescu and Fuller (2014) given by equation (15). More than one value is computed if in \emph{thetaFun} more than one population characteristic is defined.
    }
  \item{estMSE_db_telesc_WDZ}{value/s of the double bootstrap MSE estimator computed as the mean of the sums of the following elements: squared first-level boostraped error, squared first-level boostrap error for the next iteration and the opposite of second-level boostraped error based on 1 iteration (but negative sums are replaced by squared first-level boostraped error). More than one value is computed if in \emph{thetaFun} more than one population characteristic is defined.
  }
\item{estMSE_db_telesc_EF}{value/s of the telescoping double bootstrap MSE estimator proposed by Erciulescu and Fuller (2014) given by equation (15) with correction (17), where the bound for the correction is declared as \emph{q}. More than one value is computed if in \emph{thetaFun} more than one population characteristic is defined.
}

\item{estQAPE_param=estQAPE}{value/s of parametric bootstrap estimator of QAPE (Quantile of Absolute Prediction Errors) given by a quantile of absolute parametric bootstrap errors. Number of rows is equal the number of orders of quantiles to be considered (declared in \emph{p}), number of columns is equal the number of predicted characteristics (declared in  in \emph{thetaFun}).
}
\item{estQAPE_db_B2}{value/s of double-bootstrap estimator of QAPE (Quantile of Absolute Prediction Errors) given by a quantile of square roots of squared first-level boostraped errors, each corrected by the mean of squared second-level boostraped errors based on B2 iterations (where correction is made only if their difference is non-negative). Number of rows is equal the number of orders of quantiles to be considered (declared in \emph{p}), number of columns is equal the number of predicted characteristics (declared in  in \emph{thetaFun}).
}
\item{estQAPE_db_1}{value/s of double-bootstrap estimator of QAPE (Quantile of Absolute Prediction Errors) given by a quantile of square roots of squared first-level boostraped errors, each corrected by the squared second-level boostraped error based on 1 iteration (where correction is made only if their difference is non-negative). Number of rows is equal the number of orders of quantiles to be considered (declared in \emph{p}), number of columns is equal the number of predicted characteristics (declared in  in \emph{thetaFun}).

}
\item{estQAPE_db_telesc}{value/s of double-bootstrap estimator of QAPE (Quantile of Absolute Prediction Errors) given by a quantile of square roots of the sums of the following elements: squared first-level boostraped error, squared first-level boostrap error for the next iteration and the opposite of second-level boostraped error based on 1 iteration (but negative sums are replaced by squared first-level boostraped error). Number of rows is equal the number of orders of quantiles to be considered (declared in \emph{p}), number of columns is equal the number of predicted characteristics (declared in  in \emph{thetaFun}).
}

}

\references{
1. Erciulescu, A. L. and Fuller, W. A. (2013) Parametric Bootstrap Procedures for Small Area Prediction Variance. JSM 2014 - Survey Research Methods Section, pp. 3307-3318.\cr

2. Hall, P. and Maiti, T. (2006) On Parametric Bootstrap Methods for Small Area Prediction. Journal of the Royal Statistical Society. Series B, 68(2), 221-238.\cr
}
\author{
Alicja Wolny-Dominiak, Tomasz Zadlo}

\examples{

data(invData) 
#data from one period are considered: 
invData2018 <- invData[invData$year == 2018,] 
attach(invData2018)

N <- nrow(invData2018) #population size

con <- rep(1,N) 
con[c(379,380)]<-0 # last two population elements are not observed 

YS <- (investments[con==1]) # log-transformed values
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

predictor <- plugInLMM(YS, fixed.part, random.part, reg, con, weights, backTrans, thetaFun)
predictor$thetaP


# Estimation of prediction accuracy
# (q=0.77 is assumed below as in Erciulescu and Fuller (2014) eq. (17))
doubleBoot(predictor, 3, 3, c(0.5,0.9), 0.77) 
}
