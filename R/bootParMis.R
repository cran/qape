bootParMis <-
function (predictorLMM, predictorLMMmis, B, p) 
{
  if (inherits(predictorLMM, "plugInLMM") == F | inherits(predictorLMMmis, "plugInLMM") == F) {
	stop("wrong predictor/predictors")
  }
  
  if (!all.equal(predictorLMM$YS, predictorLMMmis$YS)){ 
    stop("'YS' object must be the same in both predictors")
  }
  
  if (!all.equal(predictorLMM$reg, predictorLMMmis$reg)){ 
    stop("'reg' object must be the same in both predictors")
  }
  
  if (!all.equal(predictorLMM$con, predictorLMMmis$con)){ 
    stop("'con' object must be the same in both predictors")
  }
  
  if (!all.equal(predictorLMM$thetaFun, predictorLMMmis$thetaFun)){ 
    stop("'thetaFun' object must be the same in both predictors")
  }
  
   if (B < 1) {
    stop("B must be > 0")
  }
  N <- nrow(predictorLMM$reg)
  R <- diag(predictorLMM$sigma2R/predictorLMM$weights, nrow = N, 
            ncol = N)
  namesreg <- gsub("\\.", "_", names(predictorLMM$reg))
  listVarCorr <- VarCorr(predictorLMM$mEst)
  namesvar <- sub("\\..*", "", names(listVarCorr))
  bmlist <- list()
  VarCorrNames <- attributes(listVarCorr)$names
  correct_order <- match(predictorLMM$ZBlockNames, VarCorrNames)
  for (i in 1:length(correct_order)) {
    k <- length(unique(predictorLMM$reg[, namesvar[correct_order[i]]]))
    bmlist[[i]] <- bdiag(replicate(k, as.matrix(as.data.frame(VarCorr(predictorLMM$mEst)[[correct_order[i]]])), 
                                   simplify = F))
  }
 
Gall <- bdiag(bmlist)
   if (is.positive.definite(as.matrix(Gall)) == FALSE) {
     positiveDefiniteEstG <- FALSE
     Ysim <- as.vector(predictorLMM$Xbeta) + t(chol(as.matrix(R))) %*% matrix(rnorm(N * B), nrow = N)
     cat(paste("non-positive definite estimated covariance matrix of random effects - y is generated based on a model without random effects", "\n"))
   } else 
   {positiveDefiniteEstG <- TRUE
   Ysim <- as.vector(predictorLMM$Xbeta) + 
     predictorLMM$Z %*% (t(chol(as.matrix(Gall))) %*% matrix(rnorm(ncol(predictorLMM$Z) * B), nrow = ncol(predictorLMM$Z))) + 
     t(chol(as.matrix(R))) %*% matrix(rnorm(N * B), nrow = N)}
   
   
   
   
   
  YsimS <- matrix(Ysim[predictorLMM$con == 1, ], ncol = B)

    if (is.null(predictorLMM$backTrans)) {
      predictorLMM$backTrans <- function(x) x
    }
    thetaSim <- sapply(1:B, function(i) as.numeric(predictorLMM$thetaFun(predictorLMM$backTrans(Ysim[, 
                                                                                               i]))))
    predictorLMMSim <- sapply(1:B, function(i) {
      thetaPplugin <- as.numeric(plugInLMM(YsimS[, i], 
                                           predictorLMM$fixed.part, predictorLMM$random.part, 
                                           predictorLMM$reg, predictorLMM$con, predictorLMM$weights, 
                                           predictorLMM$backTrans, predictorLMM$thetaFun)$thetaP)
      return(thetaPplugin)
    })
    
    predictorLMMmisSim <- sapply(1:B, function(i) {
      thetaPpluginmis <- as.numeric(plugInLMM(YsimS[, i], 
                                              predictorLMMmis$fixed.part, predictorLMMmis$random.part, 
                                              predictorLMMmis$reg, predictorLMMmis$con, predictorLMMmis$weights, 
                                              predictorLMMmis$backTrans, predictorLMMmis$thetaFun)$thetaP)
      return(thetaPpluginmis)
    })
    
    
  quantileNaN <- function (x, probs) {
      if (sum(is.nan(x)) > 0) rep(NaN,length(probs)) else {quantile(x, probs)}}  
    
  errorLMM <- matrix((predictorLMMSim - thetaSim), ncol = B)
  errorLMMmis <- matrix((predictorLMMmisSim - thetaSim), ncol = B)
  return(list(estQAPElmm = sapply(1:nrow(errorLMM), function(i) quantileNaN(abs(errorLMM[i, 
  ]), probs = p)), estRMSElmm = sapply(1:nrow(errorLMM), function(i) sqrt((sum(errorLMM[i, 
  ]^2))/length(errorLMM[i, ]))),
  estQAPElmmMis = sapply(1:nrow(errorLMMmis), function(i) quantileNaN(abs(errorLMMmis[i, 
  ]), probs = p)), estRMSElmmMis = sapply(1:nrow(errorLMMmis), function(i) sqrt((sum(errorLMMmis[i, 
  ]^2))/length(errorLMMmis[i, ]))),
  predictorLMMSim = predictorLMMSim, 
  predictorLMMmisSim = predictorLMMmisSim, 
  thetaSim = thetaSim, Ysim = Ysim, errorLMM = errorLMM, errorLMMmis = errorLMMmis, positiveDefiniteEstG = positiveDefiniteEstG))
}
