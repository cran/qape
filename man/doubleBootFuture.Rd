\name{doubleBootFuture}
\alias{doubleBootFuture}

\title{Double bootstrap estimators of prediction accuracy - parallel computing
}

\description{The function computes values of double bootstrap estimators of the MSE and the QAPE prediction accuracy measures using parallel computing.

}
\usage{
doubleBootFuture(predictor, B1, B2, p, q)
}

\arguments{
  \item{predictor}{one of objects: EBLUP, ebpLMMne or plugInLMM.
}
  \item{B1}{number of first-level bootstrap iterations.
}
  \item{B2}{number of second-level bootstrap iterations.
}
  \item{p}{orders of quantiles in the QAPE.
}
  \item{q}{estimator bounds assumed for \emph{estMSE_db_1_EF} and \emph{estMSE_db_telesc_EF} (which are corrected versions of \emph{estMSE_db_1} and \emph{estMSE_db_telesc}, respectively).
}
}

\details{Double-bootstrap method considered by Hall and Maiti (2006) and  Erciulescu and Fuller (2013) is used. Vectors of random effects and random components are generated from the multivariate normal distribution and REML estimates of model parameters are used. Random effects are generated for all population elements even for subsets with zero sample sizes (for which random effects are not estimated). Double-bootstrap MSE estimator presented in Hall and Maiti (2006) and  Erciulescu and Fuller (2013) are taken into account. The QAPE is a quantile of absolute prediction error which means that at least p100\% of realizations of absolute prediction errors are smaller or equal to QAPE. The parallel processing is performed via the \emph{future.apply} package.  
}


\value{

  \item{estMSE_param}{value/s of the parametric bootstrap MSE estimator. More than one value is computed if in \emph{thetaFun} more than one population characteristic is defined.
    }
  \item{estMSE_db_B2}{value/s of the double bootstrap MSE estimator computed as the difference of doubled value of \emph{estMSE_param} and the second-level MSE estimator based on B2 iterations. More than one value is computed if in \emph{thetaFun} more than one population characteristic is defined.
    }
  \item{estMSE_db_B2_WDZ}{value/s of the double bootstrap MSE estimator computed as the mean of squared first-level bootstraped errors, each corrected by the mean of squared second-level bootstraped errors based on B2 iterations (where correction is made only if their difference is non-negative). More than one value is computed if in \emph{thetaFun} more than one population characteristic is defined.
    }
  \item{estMSE_db_B2_HM}{value/s of the double bootstrap MSE estimator proposed by Hall and Maiti (2006) equation (2.17). More than one value is computed if in \emph{thetaFun} more than one population characteristic is defined.
    }
  \item{estMSE_db_1}{value/s of the double bootstrap MSE estimator computed as the difference of doubled value of \emph{estMSE_param} and the second-level MSE estimator based on B2=1 iteration. More than one value is computed if in \emph{thetaFun} more than one population characteristic is defined.
    }
  
  \item{estMSE_db_1_WDZ}{value/s of the double bootstrap MSE estimator computed as the mean of squared first-level bootstraped errors, each corrected by the squared second-level bootstraped error based on 1 iteration (where correction is made only if their difference is non-negative). More than one value is computed if in \emph{thetaFun} more than one population characteristic is defined.
    }
  \item{estMSE_db_1_EF}{value/s of the double bootstrap MSE estimator proposed by Erciulescu and Fuller (2014) given by equation (13) with correction (17), where the bound for the correction is declared as \emph{q}. More than one value is computed if in \emph{thetaFun} more than one population characteristic is defined.
    }
  \item{estMSE_db_telesc}{value/s of the telescoping double bootstrap MSE estimator proposed by Erciulescu and Fuller (2014) given by equation (15). More than one value is computed if in \emph{thetaFun} more than one population characteristic is defined.
    }
  \item{estMSE_db_telesc_WDZ}{value/s of the double bootstrap MSE estimator computed as the mean of the sums of the following elements: squared first-level bootstraped error, squared first-level bootstrap error for the next iteration and the opposite of second-level bootstraped error based on 1 iteration (but negative sums are replaced by squared first-level bootstraped error). More than one value is computed if in \emph{thetaFun} more than one population characteristic is defined.
  }
\item{estMSE_db_telesc_EF}{value/s of the telescoping double bootstrap MSE estimator proposed by Erciulescu and Fuller (2014) given by equation (15) with correction (17), where the bound for the correction is declared as \emph{q}. More than one value is computed if in \emph{thetaFun} more than one population characteristic is defined.
}

\item{estQAPE_param}{value/s of parametric bootstrap estimator of QAPE (Quantile of Absolute Prediction Error) given by a quantile of absolute parametric bootstrap errors. Number of rows is equal the number of orders of quantiles to be considered (declared in \emph{p}), number of columns is equal the number of predicted characteristics (declared in  in \emph{thetaFun}).
}

\item{estQAPE_db_B2}{value/s of double-bootstrap estimator of QAPE (Quantile of Absolute Prediction Error) given by a quantile of square roots of squared first-level bootstraped errors, each corrected by the mean of squared second-level bootstraped errors based on B2 iterations (where correction is made only if their difference is non-negative). Number of rows is equal the number of orders of quantiles to be considered (declared in \emph{p}), number of columns is equal the number of predicted characteristics (declared in  in \emph{thetaFun}).
}
\item{estQAPE_db_1}{value/s of double-bootstrap estimator of QAPE (Quantile of Absolute Prediction Error) given by a quantile of square roots of squared first-level bootstraped errors, each corrected by the squared second-level bootstraped error based on 1 iteration (where correction is made only if their difference is non-negative). Number of rows is equal the number of orders of quantiles to be considered (declared in \emph{p}), number of columns is equal the number of predicted characteristics (declared in  in \emph{thetaFun}).

}
\item{estQAPE_db_telesc}{value/s of double-bootstrap estimator of QAPE (Quantile of Absolute Prediction Error) given by a quantile of square roots of the sums of the following elements: squared first-level bootstraped error, squared first-level bootstrap error for the next iteration and the opposite of second-level bootstraped error based on 1 iteration (but negative sums are replaced by squared first-level bootstraped error). Number of rows is equal to the number of orders of quantiles to be considered (declared in \emph{p}), number of columns is equal to the number of predicted characteristics (declared in  in \emph{thetaFun}).
}
\item{error1}{the matrix of first-level bootstrap errors. Number of rows is equal to the number of predicted characteristics (declared in  in \emph{thetaFun}), number of columns is equal to \emph{B1}.
}
\item{error2}{the list of matrices of second-level bootstrap errors. The length of list is equal to the number of predicted characteristics (declared in  in \emph{thetaFun}), the number of rows of each matrix is equal to \emph{B1}, the number of columns is equal to \emph{B2}.
}
\item{corSquaredError1_db_B2}{the matrix of corrected squared first-level bootstrap errors defined as doubled squared first-level bootstrap errors minus the mean of squared second-level bootstrap errors (computed for the approriate first-level bootstrap iterations). Number of rows is equal to \emph{B1}, the number of columns is equal to the number of predicted characteristics (declared in  in \emph{thetaFun}). Values can be negative.
}
\item{corSquaredError1_db_1}{the matrix of corrected squared first-level bootstrap errors defined as doubled squared first-level bootstrap errors minus the squared second-level bootstrap error (computed once for each first-level bootstrap iteration). Number of rows is equal to \emph{B1}, the number of columns is equal to the number of predicted characteristics (declared in  in \emph{thetaFun}). Values can be negative.
}

\item{corSquaredError1_db_telesc}{the matrix of corrected squared first-level bootstrap errors defined by elements from which the average given by equation (15) in Erciulescu and Fuller (2014) is counted. Number of rows is equal to \emph{B1}, the number of columns is equal to the number of predicted characteristics (declared in  in \emph{thetaFun}). Values can be negative.
}

\item{corSquaredError1_db_B2_WDZ}{the matrix of squared first-level bootstraped errors, each corrected by the mean of squared second-level bootstraped errors based on B2 iterations (where correction is made only if their difference is non-negative). Number of rows is equal to \emph{B1}, the number of columns is equal to the number of predicted characteristics (declared in  in \emph{thetaFun}). Values are non-negative.}

\item{corSquaredError1_db_1_WDZ}{the matrix of squared first-level bootstraped errors, each corrected by the squared second-level bootstraped error based on 1 iteration (where correction is made only if their difference is non-negative). Number of rows is equal to \emph{B1}, the number of columns is equal to the number of predicted characteristics (declared in  in \emph{thetaFun}). Values are non-negative.}

\item{corSquaredError1_db_telesc_WDZ}{the matrix of corrected squared first-level bootstrap errors defined by sums of the following elements: squared first-level bootstraped error, squared first-level bootstrap error for the next iteration and the opposite of second-level bootstraped error based on 1 iteration (but negative sums are replaced by squared first-level bootstraped error).Number of rows is equal to \emph{B1}, the number of columns is equal to the number of predicted characteristics (declared in  in \emph{thetaFun}). Values are non-negative.}


\item{positiveDefiniteEstGlev1}{logical indicating if the estimated covariance matrix of random effects, used to generate bootstrap realizations of the dependent variable at the first level of the double bootstrap, is positive definite.} 

\item{positiveDefiniteEstGlev2}{number of cases ouf of B1 with positive definite estimated covariance matrix of random effects used to generate bootstrap realizations of the dependent variable at the second level of the double bootstrap.}
}

\references{
1. Erciulescu, A. L. and Fuller, W. A. (2013) Parametric Bootstrap Procedures for Small Area Prediction Variance. JSM 2014 - Survey Research Methods Section, 3307-3318.\cr

2. Hall, P. and Maiti, T. (2006) On Parametric Bootstrap Methods for Small Area Prediction. Journal of the Royal Statistical Society. Series B, 68(2), 221-238.\cr
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

### Characteristics to be predicted:
# values of the variable for last two population elements  
thetaFun <- function(x) {x[c(379,380)]}
set.seed(123456)

predictor <- plugInLMM(YS, fixed.part, random.part, reg, con, weights, backTrans, thetaFun)
predictor$thetaP


### Estimation of prediction accuracy
# in the first column
# for the predictor of the value of the variable for population element no. 379,
# in the second column
# for the predictor of the value of the variable for population element no. 380:
doubleBootFuture(predictor, 3, 3, c(0.5,0.9), 0.77) 
#q=0.77 assumed as in Erciulescu and FUller (2014) eq. (17)

detach(invData2018)
}

