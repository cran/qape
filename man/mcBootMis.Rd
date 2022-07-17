\name{mcBootMis}
\alias{mcBootMis}
\title{Monte Carlo simulation study of accuracy of estimators of accuracy measures
}
\description{The function computes in the Monte Carlo simulation study values of accuracy measures of estimators of accuracy measures of two predictors under the model defined by the first of them.}

\usage{
mcBootMis(Ypop, predictorLMM, predictorLMMmis, K, B1, B2, p, q)
}

\arguments{
  \item{Ypop}{population values of the dependent variable.
}
  \item{predictorLMM}{plugInLMM object, the predictor used to define the model assumed in the simulation study.
}
  \item{predictorLMMmis}{plugInLMM object, the second predictor, the properties of which are assessed under the misspecified model used in \emph{predictorLMM}.
}
 \item{K}{the number of Monte Carlo iterations.
}
  \item{B1}{the number of first-level bootstrap iterations.
}
  \item{B2}{the number of second-level bootstrap iterations.
}
  \item{p}{orders of quantiles in the QAPE.
}
  \item{q}{estimator bounds assumed for \emph{estMSE_db_1_EF} and \emph{estMSE_db_telesc_EF} (which are corrected versions of \emph{estMSE_db_1} and \emph{estMSE_db_telesc}, respectively).
}
}

\details{In the model-based simulation study population values of the dependent variable are generated based on the (possibly transformed) Linear Mixed Model used in \emph{predictorLMM} and the accuracy of predictors \emph{predictorLMM} and \emph{predictorLMMmis} is assessed. What is more, the the accuracy of parametric, residual and double bootstrap estimators of accuracy measures is studied under the model used in \emph{predictorLMM}. Values of some MSE estimators can be negative, the number of negative values of MSE estimators obtained in the simulation study are presented in objects \emph{neg_estMSE_LMM} and \emph{neg_estMSE_LMMmis}. Hence, some RMSE estimators computed as square roots of MSE estimators can produce NaNs - see warnings.}     


\value{
\item{QAPElmm}{value/s of the QAPE of \emph{predictorLMM} assessed in the Monte Carlo study - the number of rows is equal to the number of orders of quantiles to be considered (declared in \emph{p}), the number of columns is equal to the number of predicted characteristics (declared in \emph{thetaFun}).}

\item{RMSElmm}{value/s of the RMSE of \emph{predictorLMM} assessed in the Monte Carlo study (more than one value is computed if in \emph{thetaFun} more than one population characteristic is defined).}

\item{rRMSElmm}{value/s of the rRMSE (in percentages) of \emph{predictorLMM} assessed in the Monte Carlo study (more than one value is computed if in \emph{thetaFun} more than one population characteristic is defined).}

\item{rBlmm}{value/s of the relative bias (in percentages) of \emph{predictorLMM} assessed in the Monte Carlo study (more than one value is computed if in \emph{thetaFun} more than one population characteristic is defined).}

\item{QAPElmmMis}{value/s of the QAPE of \emph{predictorLMM2} assessed in the Monte Carlo study - the number of rows is equal to the number of orders of quantiles to be considered (declared in \emph{p}), the number of columns is equal to the number of predicted characteristics (declared in \emph{thetaFun}).}

\item{RMSElmmMis}{value/s of the RMSE of \emph{predictorLMM2} assessed in the Monte Carlo study (more than one value is computed if in \emph{thetaFun} more than one population characteristic is defined).}

\item{rRMSElmmMis}{value/s of the rRMSE (in percentages) of \emph{predictorLMM2} assessed in the Monte Carlo study (more than one value is computed if in \emph{thetaFun} more than one population characteristic is defined).}

\item{rBlmmMis}{value/s of the relative bias (in percentages) of \emph{predictorLMMmis} assessed in the Monte Carlo study (more than one value is computed if in \emph{thetaFun} more than one population characteristic is defined).}

\item{rB.estRMSE_rbF_LMM}{relative bias (in percentages) of estimated value/s of RMSE of \emph{predictorLMM} without correction to avoid the problem of underdispersion of residual bootstrap distributions  (more than one value is computed if in \emph{thetaFun} more than one population characteristic is defined).}

\item{rRMSE.estRMSE_rbF_LMM}{relative RMSE (in percentages) of estimated value/s of RMSE of \emph{predictorLMM} without correction to avoid the problem of underdispersion of residual bootstrap distributions  (more than one value is computed if in \emph{thetaFun} more than one population characteristic is defined).}

\item{rB.estRMSE_rbF_LMMmis}{relative bias (in percentages) of estimated value/s of RMSE of \emph{predictorLMMmis} without correction to avoid the problem of underdispersion of residual bootstrap distributions  (more than one value is computed if in \emph{thetaFun} more than one population characteristic is defined).}

\item{rRMSE.estRMSE_rbF_LMMmis}{relative RMSE (in percentages) of estimated value/s of RMSE of \emph{predictorLMMmis} without correction to avoid the problem of underdispersion of residual bootstrap distributions  (more than one value is computed if in \emph{thetaFun} more than one population characteristic is defined).}


\item{rB.estMSE_rbF_LMM}{relative bias (in percentages) of estimated value/s of MSE of \emph{predictorLMM} without correction to avoid the problem of underdispersion of residual bootstrap distributions  (more than one value is computed if in \emph{thetaFun} more than one population characteristic is defined).}

\item{rRMSE.estMSE_rbF_LMM}{relative RMSE (in percentages) of estimated value/s of MSE of \emph{predictorLMM} without correction to avoid the problem of underdispersion of residual bootstrap distributions  (more than one value is computed if in \emph{thetaFun} more than one population characteristic is defined).}

\item{rB.estMSE_rbF_LMMmis}{relative bias (in percentages) of estimated value/s of MSE of \emph{predictorLMMmis} without correction to avoid the problem of underdispersion of residual bootstrap distributions  (more than one value is computed if in \emph{thetaFun} more than one population characteristic is defined).}

\item{rRMSE.estMSE_rbF_LMMmis}{relative RMSE (in percentages) of estimated value/s of MSE of \emph{predictorLMMmis} without correction to avoid the problem of underdispersion of residual bootstrap distributions  (more than one value is computed if in \emph{thetaFun} more than one population characteristic is defined).}

\item{rB.estQAPE_rbF_LMM}{relative bias (in percentages) of estimated value/s of QAPE of \emph{predictorLMM} without correction to avoid the problem of underdispersion of residual bootstrap distributions, the number of rows is equal to the number of orders of quantiles to be considered (declared in \emph{p}), the number of columns is equal to the number of predicted characteristics (declared in \emph{thetaFun}).}

\item{rRMSE.estQAPE_rbF_LMM}{relative RMSE (in percentages) of estimated value/s of QAPE of \emph{predictorLMM} without correction to avoid the problem of underdispersion of residual bootstrap distributions, the number of rows is equal to the number of orders of quantiles to be considered (declared in \emph{p}), the number of columns is equal to the number of predicted characteristics (declared in \emph{thetaFun}).}

\item{rB.estQAPE_rbF_LMMmis}{relative bias (in percentages) of estimated value/s of QAPE of \emph{predictorLMMmis} without correction to avoid the problem of underdispersion of residual bootstrap distributions, the number of rows is equal to the number of orders of quantiles to be considered (declared in \emph{p}), the number of columns is equal to the number of predicted characteristics (declared in \emph{thetaFun}).}

\item{rRMSE.estQAPE_rbF_LMMmis}{relative RMSE (in percentages) of estimated value/s of QAPE of \emph{predictorLMMmis} without correction to avoid the problem of underdispersion of residual bootstrap distributions, the number of columns is equal to the number of predicted characteristics (declared in \emph{thetaFun}) .}

\item{rB.estRMSE_rbT_LMM}{relative bias (in percentages) of estimated value/s of RMSE of \emph{predictorLMM} with correction to avoid the problem of underdispersion of residual bootstrap distributions  (more than one value is computed if in \emph{thetaFun} more than one population characteristic is defined).}

\item{rRMSE.estRMSE_rbT_LMM}{relative RMSE (in percentages) of estimated value/s of RMSE of \emph{predictorLMM} with correction to avoid the problem of underdispersion of residual bootstrap distributions  (more than one value is computed if in \emph{thetaFun} more than one population characteristic is defined).}

\item{rB.estRMSE_rbT_LMMmis}{relative bias (in percentages) of estimated value/s of RMSE of \emph{predictorLMMmis} with correction to avoid the problem of underdispersion of residual bootstrap distributions  (more than one value is computed if in \emph{thetaFun} more than one population characteristic is defined).}

\item{rRMSE.estRMSE_rbT_LMMmis}{relative RMSE (in percentages) of estimated value/s of RMSE of \emph{predictorLMMmis} with correction to avoid the problem of underdispersion of residual bootstrap distributions  (more than one value is computed if in \emph{thetaFun} more than one population characteristic is defined).}

\item{rB.estMSE_rbT_LMM}{relative bias (in percentages) of estimated value/s of MSE of \emph{predictorLMM} with correction to avoid the problem of underdispersion of residual bootstrap distributions  (more than one value is computed if in \emph{thetaFun} more than one population characteristic is defined).}

\item{rRMSE.estMSE_rbT_LMM}{relative RMSE (in percentages) of estimated value/s of MSE of \emph{predictorLMM} with correction to avoid the problem of underdispersion of residual bootstrap distributions  (more than one value is computed if in \emph{thetaFun} more than one population characteristic is defined).}

\item{rB.estMSE_rbT_LMMmis}{relative bias (in percentages) of estimated value/s of MSE of \emph{predictorLMMmis} with correction to avoid the problem of underdispersion of residual bootstrap distributions  (more than one value is computed if in \emph{thetaFun} more than one population characteristic is defined).}

\item{rRMSE.estMSE_rbT_LMMmis}{relative RMSE (in percentages) of estimated value/s of MSE of \emph{predictorLMMmis} with correction to avoid the problem of underdispersion of residual bootstrap distributions  (more than one value is computed if in \emph{thetaFun} more than one population characteristic is defined).}

\item{rB.estQAPE_rbT_LMM}{relative bias (in percentages) of estimated value/s of QAPE of \emph{predictorLMM} with correction to avoid the problem of underdispersion of residual bootstrap distributions  (more than one value is computed if in \emph{thetaFun} more than one population characteristic is defined).}

\item{rRMSE.estQAPE_rbT_LMM}{relative RMSE (in percentages) of estimated value/s of QAPE of \emph{predictorLMM} with correction to avoid the problem of underdispersion of residual bootstrap distributions  (more than one value is computed if in \emph{thetaFun} more than one population characteristic is defined).}

\item{rB.estQAPE_rbT_LMMmis}{relative bias (in percentages) of estimated value/s of QAPE of \emph{predictorLMMmis} with correction to avoid the problem of underdispersion of residual bootstrap distributions  (more than one value is computed if in \emph{thetaFun} more than one population characteristic is defined).}

\item{rRMSE.estQAPE_rbT_LMMmis}{relative RMSE (in percentages) of estimated value/s of QAPE of \emph{predictorLMMmis} with correction to avoid the problem of underdispersion of residual bootstrap distributions  (more than one value is computed if in \emph{thetaFun} more than one population characteristic is defined).}


\item{neg_estMSE_LMM}{the number of negative values of MSE estimators of \emph{predictorLMM} obtained in the simulaton study out of K iterations, the number of rows is equal to 10 - the number of considered parametric and double bootstrap MSE estimators, the number of columns is equal to the number of predicted characteristics (declared in \emph{thetaFun}).}

\item{neg_estMSE_LMMmis}{the number of negative values of MSE estimators of \emph{predictorLMMmis} obtained in the simulaton study out of K iterations, the number of rows is equal to 10 - the number of considered parametric and double bootstrap MSE estimators, the number of columns is equal to the number of predicted characteristics (declared in \emph{thetaFun}).}

\item{rB.estMSE_param_LMMmis}{relative bias (in percentages) of \emph{estMSE_param} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.}

\item{rRMSE.estMSE_param_LMMmis}{relative RMSE (in percentages) of \emph{estMSE_param} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.}

\item{rB.estMSE_db_B2_LMMmis}{relative bias (in percentages) of \emph{estMSE_db_B2} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.}

\item{rRMSE.estMSE_db_B2_LMMmis}{relative RMSE (in percentages) of \emph{estMSE_db_B2} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.}

\item{rB.estMSE_db_B2_WDZ_LMMmis}{relative bias (in percentages) of \emph{estMSE_db_B2_WDZ} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.}

\item{rRMSE.estMSE_db_B2_WDZ_LMMmis}{relative RMSE (in percentages) of \emph{estMSE_db_B2_WDZ} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.}

\item{rB.estMSE_db_B2_HM_LMMmis}{relative bias (in percentages) of \emph{estMSE_db_B2_HM} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.}

\item{rRMSE.estMSE_db_B2_HM_LMMmis}{relative RMSE (in percentages) of \emph{estMSE_db_B2_HM} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.}

\item{rB.estMSE_db_1_LMMmis}{relative bias (in percentages) of \emph{estMSE_db_1} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.}

\item{rRMSE.estMSE_db_1_LMMmis}{relative RMSE (in percentages) of \emph{estMSE_db_1} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.}

\item{rB.estMSE_db_1_WDZ_LMMmis}{relative bias (in percentages) of \emph{estMSE_db_1_WDZ} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.}

\item{rRMSE.estMSE_db_1_WDZ_LMMmis}{relative RMSE (in percentages) of \emph{estMSE_db_1_WDZ} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.}

\item{rB.estMSE_db_1_EF_LMMmis}{relative bias (in percentages) of \emph{estMSE_db_1_EF} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.}

\item{rRMSE.estMSE_db_1_EF_LMMmis}{relative RMSE (in percentages) of \emph{estMSE_db_1_EF} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.}

\item{rB.estMSE_db_telesc_LMMmis}{relative bias (in percentages) of \emph{estMSE_db_telesc} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.}

\item{rRMSE.estMSE_db_telesc_LMMmis}{relative RMSE (in percentages) of \emph{estMSE_db_telesc} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.}

\item{rB.estMSE_db_telesc_WDZ_LMMmis}{relative bias (in percentages) of \emph{estMSE_db_telesc_WDZ} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.}

\item{rRMSE.estMSE_db_telesc_WDZ_LMMmis}{relative RMSE (in percentages) of \emph{estMSE_db_telesc_WDZ} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.}

\item{rB.estMSE_db_telesc_EF_LMMmis}{relative bias (in percentages) of \emph{estMSE_db_telesc_EF} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.}

\item{rRMSE.estMSE_db_telesc_EF_LMMmis}{relative RMSE (in percentages) of \emph{estMSE_db_telesc_EF} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.}

\item{rB.estRMSE_param_LMMmis}{relative bias (in percentages) of \emph{estRMSE_param} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.}

\item{rRMSE.estRMSE_param_LMMmis}{relative RMSE (in percentages) of \emph{estRMSE_param} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.}

\item{rB.estRMSE_db_B2_LMMmis}{relative bias (in percentages) of \emph{estRMSE_db_B2} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.}

\item{rRMSE.estRMSE_db_B2_LMMmis}{relative RMSE (in percentages) of \emph{estRMSE_db_B2} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.}

\item{rB.estRMSE_db_B2_WDZ_LMMmis}{relative bias (in percentages) of \emph{estRMSE_db_B2_WDZ} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.}

\item{rRMSE.estRMSE_db_B2_WDZ_LMMmis}{relative RMSE (in percentages) of \emph{estRMSE_db_B2_WDZ} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.}

\item{rB.estRMSE_db_B2_HM_LMMmis}{relative bias (in percentages) of \emph{estRMSE_db_B2_HM} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.}

\item{rRMSE.estRMSE_db_B2_HM_LMMmis}{relative RMSE (in percentages) of \emph{estRMSE_db_B2_HM} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.}

\item{rB.estRMSE_db_1_LMMmis}{relative bias (in percentages) of \emph{estRMSE_db_1} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.}

\item{rRMSE.estRMSE_db_1_LMMmis}{relative RMSE (in percentages) of \emph{estRMSE_db_1} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.}

\item{rB.estRMSE_db_1_WDZ_LMMmis}{relative bias (in percentages) of \emph{estRMSE_db_1_WDZ} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.}

\item{rRMSE.estRMSE_db_1_WDZ_LMMmis}{relative RMSE (in percentages) of \emph{estRMSE_db_1_WDZ} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.}

\item{rB.estRMSE_db_1_EF_LMMmis}{relative bias (in percentages) of \emph{estRMSE_db_1_EF} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.}

\item{rRMSE.estRMSE_db_1_EF_LMMmis}{relative RMSE (in percentages) of \emph{estRMSE_db_1_EF} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.}

\item{rB.estRMSE_db_telesc_LMMmis}{relative bias (in percentages) of \emph{estRMSE_db_telesc} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.}

\item{rRMSE.estRMSE_db_telesc_LMMmis}{relative RMSE (in percentages) of \emph{estRMSE_db_telesc} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.}

\item{rB.estRMSE_db_telesc_WDZ_LMMmis}{relative bias (in percentages) of \emph{estRMSE_db_telesc_WDZ} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.}

\item{rRMSE.estRMSE_db_telesc_WDZ_LMMmis}{relative RMSE (in percentages) of \emph{estRMSE_db_telesc_WDZ} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.}

\item{rB.estRMSE_db_telesc_EF_LMMmis}{relative bias (in percentages) of \emph{estRMSE_db_telesc_EF} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.}

\item{rRMSE.estRMSE_db_telesc_EF_LMMmis}{relative RMSE (in percentages) of \emph{estRMSE_db_telesc_EF} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.}

\item{rB.estQAPE_param_LMMmis}{relative bias (in percentages) of \emph{estQAPE_param} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.}

\item{rRMSE.estQAPE_param_LMMmis}{relative RMSE (in percentages) of \emph{estQAPE_param} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.}

\item{rB.estQAPE_db_B2_LMMmis}{relative bias (in percentages) of \emph{estQAPE_db_B2} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.} 

\item{rRMSE.estQAPE_db_B2_LMMmis}{relative RMSE (in percentages) of \emph{estQAPE_db_B2} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.} 

\item{rB.estQAPE_db_1_LMMmis}{relative bias (in percentages) of \emph{estQAPE_db_1} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.}  

\item{rRMSE.estQAPE_db_1_LMMmis}{relative RMSE (in percentages) of \emph{estQAPE_db_1} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.}  

\item{rB.estQAPE_db_telesc_LMMmis}{relative bias (in percentages) of \emph{estQAPE_db_telesc} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.}  

\item{rRMSE.estQAPE_db_telesc_LMMmis}{relative RMSE (in percentages) of \emph{estQAPE_db_telesc} estimator (see \emph{doubleBoot} function) of \emph{predictorLMMmis}.}   

\item{rB.estMSE_param_LMM}{relative bias (in percentages) of \emph{estMSE_param} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.}

\item{rRMSE.estMSE_param_LMM}{relative RMSE (in percentages) of \emph{estMSE_param} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.}

\item{rB.estMSE_db_B2_LMM}{relative bias (in percentages) of \emph{estMSE_db_B2} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.}

\item{rRMSE.estMSE_db_B2_LMM}{relative RMSE (in percentages) of \emph{estMSE_db_B2} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.}

\item{rB.estMSE_db_B2_WDZ_LMM}{relative bias (in percentages) of \emph{estMSE_db_B2_WDZ} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.}

\item{rRMSE.estMSE_db_B2_WDZ_LMM}{relative RMSE (in percentages) of \emph{estMSE_db_B2_WDZ} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.}

\item{rB.estMSE_db_B2_HM_LMM}{relative bias (in percentages) of \emph{estMSE_db_B2_HM} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.}

\item{rRMSE.estMSE_db_B2_HM_LMM}{relative RMSE (in percentages) of \emph{estMSE_db_B2_HM} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.}

\item{rB.estMSE_db_1_LMM}{relative bias (in percentages) of \emph{estMSE_db_1} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.}

\item{rRMSE.estMSE_db_1_LMM}{relative RMSE (in percentages) of \emph{estMSE_db_1} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.}

\item{rB.estMSE_db_1_WDZ_LMM}{relative bias (in percentages) of \emph{estMSE_db_1_WDZ} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.}

\item{rRMSE.estMSE_db_1_WDZ_LMM}{relative RMSE (in percentages) of \emph{estMSE_db_1_WDZ} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.}

\item{rB.estMSE_db_1_EF_LMM}{relative bias (in percentages) of \emph{estMSE_db_1_EF} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.}

\item{rRMSE.estMSE_db_1_EF_LMM}{relative RMSE (in percentages) of \emph{estMSE_db_1_EF} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.}

\item{rB.estMSE_db_telesc_LMM}{relative bias (in percentages) of \emph{estMSE_db_telesc} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.}

\item{rRMSE.estMSE_db_telesc_LMM}{relative RMSE (in percentages) of \emph{estMSE_db_telesc} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.}

\item{rB.estMSE_db_telesc_WDZ_LMM}{relative bias (in percentages) of \emph{estMSE_db_telesc_WDZ} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.}

\item{rRMSE.estMSE_db_telesc_WDZ_LMM}{relative RMSE (in percentages) of \emph{estMSE_db_telesc_WDZ} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.}

\item{rB.estMSE_db_telesc_EF_LMM}{relative bias (in percentages) of \emph{estMSE_db_telesc_EF} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.}

\item{rRMSE.estMSE_db_telesc_EF_LMM}{relative RMSE (in percentages) of \emph{estMSE_db_telesc_EF} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.}

\item{rB.estRMSE_param_LMM}{relative bias (in percentages) of \emph{estRMSE_param} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.}

\item{rRMSE.estRMSE_param_LMM}{relative RMSE (in percentages) of \emph{estRMSE_param} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.}

\item{rB.estRMSE_db_B2_LMM}{relative bias (in percentages) of \emph{estRMSE_db_B2} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.}

\item{rRMSE.estRMSE_db_B2_LMM}{relative RMSE (in percentages) of \emph{estRMSE_db_B2} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.}

\item{rB.estRMSE_db_B2_WDZ_LMM}{relative bias (in percentages) of \emph{estRMSE_db_B2_WDZ} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.}

\item{rRMSE.estRMSE_db_B2_WDZ_LMM}{relative RMSE (in percentages) of \emph{estRMSE_db_B2_WDZ} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.}

\item{rB.estRMSE_db_B2_HM_LMM}{relative bias (in percentages) of \emph{estRMSE_db_B2_HM} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.}

\item{rRMSE.estRMSE_db_B2_HM_LMM}{relative RMSE (in percentages) of \emph{estRMSE_db_B2_HM} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.}

\item{rB.estRMSE_db_1_LMM}{relative bias (in percentages) of \emph{estRMSE_db_1} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.}

\item{rRMSE.estRMSE_db_1_LMM}{relative RMSE (in percentages) of \emph{estRMSE_db_1} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.}

\item{rB.estRMSE_db_1_WDZ_LMM}{relative bias (in percentages) of \emph{estRMSE_db_1_WDZ} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.}

\item{rRMSE.estRMSE_db_1_WDZ_LMM}{relative RMSE (in percentages) of \emph{estRMSE_db_1_WDZ} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.}

\item{rB.estRMSE_db_1_EF_LMM}{relative bias (in percentages) of \emph{estRMSE_db_1_EF} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.}

\item{rRMSE.estRMSE_db_1_EF_LMM}{relative RMSE (in percentages) of \emph{estRMSE_db_1_EF} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.}

\item{rB.estRMSE_db_telesc_LMM}{relative bias (in percentages) of \emph{estRMSE_db_telesc} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.}

\item{rRMSE.estRMSE_db_telesc_LMM}{relative RMSE (in percentages) of \emph{estRMSE_db_telesc} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.}

\item{rB.estRMSE_db_telesc_WDZ_LMM}{relative bias (in percentages) of \emph{estRMSE_db_telesc_WDZ} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.}

\item{rRMSE.estRMSE_db_telesc_WDZ_LMM}{relative RMSE (in percentages) of \emph{estRMSE_db_telesc_WDZ} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.}

\item{rB.estRMSE_db_telesc_EF_LMM}{relative bias (in percentages) of \emph{estRMSE_db_telesc_EF} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.}

\item{rRMSE.estRMSE_db_telesc_EF_LMM}{relative RMSE (in percentages) of \emph{estRMSE_db_telesc_EF} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.}

\item{rB.estQAPE_param_LMM}{relative bias (in percentages) of \emph{estQAPE_param} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.}

\item{rRMSE.estQAPE_param_LMM}{relative RMSE (in percentages) of \emph{estQAPE_param} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.}

\item{rB.estQAPE_db_B2_LMM}{relative bias (in percentages) of \emph{estQAPE_db_B2} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.} 

\item{rRMSE.estQAPE_db_B2_LMM}{relative RMSE (in percentages) of \emph{estQAPE_db_B2} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.} 

\item{rB.estQAPE_db_1_LMM}{relative bias (in percentages) of \emph{estQAPE_db_1} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.}  

\item{rRMSE.estQAPE_db_1_LMM}{relative RMSE (in percentages) of \emph{estQAPE_db_1} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.}  

\item{rB.estQAPE_db_telesc_LMM}{relative bias (in percentages) of \emph{estQAPE_db_telesc} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.}  

\item{rRMSE.estQAPE_db_telesc_LMM}{relative RMSE (in percentages) of \emph{estQAPE_db_telesc} estimator (see \emph{doubleBoot} function) of \emph{predictorLMM}.}  

\item{MCpositiveDefiniteEstGlev1}{number of cases ouf of K with postive definite estimated covariance matrix of random effects used to generate bootstrap realizations of the dependent variable at the first level of the double bootstrap.} 

\item{MCpositiveDefiniteEstGlev2}{number of cases ouf of K*B1 with positive definite estimated covariance matrix of random effects used to generate bootstrap realizations of the dependent variable at the second level of the double bootstrap.} 
}


\references{

1. Chatterjee, S., Lahiri, P. Li, H. (2008) Parametric bootstrap approximation to the distribution of EBLUP and related prediction intervals in linear mixed models, Annals of Statistics, Vol. 36 (3), pp. 1221?1245.\cr

2. Rao, J.N.K. and Molina, I. (2015) Small Area Estimation. Second edition, John Wiley & Sons, New Jersey.\cr

3. Zadlo T. (2017), On asymmetry of prediction errors in small area estimation, Statistics in Transition, 18 (3), 413-432. \cr

}

\author{
Tomasz Zadlo, Alicja Wolny-Dominiak}

\examples{
library(lme4)


data(sleepstudy)
sleepstudy$Day.sim <- sample(sleepstudy$Days, nrow(sleepstudy)) 

Ypop <- sleepstudy$Reaction
con <- ifelse(sleepstudy$Days != 9, 1, 0)
YS <- sleepstudy$Reaction[con == 1]
fixed.part <- 'Days'
random.part <- '(1 | Subject)' 
reg <- sleepstudy[ , - which(colnames(sleepstudy)== 'Reaction')]
gamma <-  rep(1, nrow(reg))
weights <- rep(2, nrow(reg))
estMSE <- TRUE

# PLUG-IN predictor
thetaFun <- function(x) {c(sum(x), mean(x))}
plugin <- plugInLMM(YS, fixed.part, random.part, reg, 
                     con, weights, backTrans = NULL, thetaFun)

fixed.part.mis <- 'Day.sim'
plugin.mis<-plugInLMM(YS, fixed.part.mis, random.part, reg, con, weights, backTrans=NULL,thetaFun)

mcBootMis(Ypop, plugin, plugin.mis, 2, 2, 1, c(0.5, 0.9), 0.77)

}




