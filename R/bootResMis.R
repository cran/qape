bootResMis <-
function(predictorLMM, predictorLMMmis, B, p, correction) 
{
  N <- nrow(predictorLMM$X)
  YsimF <- function(...) {
    tablsrswrRe <- srswrRe(ranef(predictorLMM$mEst), predictorLMM$reg)$tablsrswrRe
    dfeS <- data.frame(nrw = 1:nrow(predictorLMM$regS), predictorLMM$regS, 
                       eS = predictorLMM$eS)
    dfesamp <- dfeS[sample(nrow(dfeS), N, replace = T), ]
    vsample <- tablsrswrRe$ranef
    names(vsample) <- tablsrswrRe$refNames
    return(predictorLMM$Xbeta + predictorLMM$Z %*% as.matrix(vsample[colnames(predictorLMM$Z)]) + 
             dfesamp$eS)
  }
  YsimT <- function(...) {
    lranefKorekta <- corrRanef(predictorLMM$mEst)
    tablwzlKorekta <- srswrRe(lranefKorekta, predictorLMM$reg)$tablsrswrRe
    eSKorekta <- corrRancomp(predictorLMM$mEst)
    dfeSKorekta <- data.frame(nrw = 1:nrow(predictorLMM$regS), 
                              predictorLMM$regS, eSKorekta)
    dfesampKorekta <- dfeSKorekta[sample(dfeSKorekta$nrw, 
                                         N, replace = T), ]
    vsampleKorekta <- tablwzlKorekta$ranef
    names(vsampleKorekta) <- tablwzlKorekta$refNames
    as.matrix(vsampleKorekta[colnames(predictorLMM$Z)])
    return(predictorLMM$Xbeta + predictorLMM$Z %*% as.matrix(vsampleKorekta[colnames(predictorLMM$Z)]) + 
             dfesampKorekta$eSKorekta)
  }
  
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
  if (correction == F) {
    Ysim <- matrix(replicate(B, YsimF(predictorLMM)), ncol = B)
  }  else {Ysim <- matrix(replicate(B, YsimT(predictorLMM)), ncol = B)}
  
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
    thetaPplugin <- as.numeric(plugInLMM(YsimS[, i], 
                                         predictorLMM$fixed.part, predictorLMM$random.part, 
                                         predictorLMM$reg, predictorLMM$con, predictorLMM$weights, 
                                         predictorLMM$backTrans, predictorLMM$thetaFun)$thetaP)
    return(thetaPplugin)
  })
  
  predictorLMMmisSim <- sapply(1:B, function(i) {
    thetaPpluginMis <- as.numeric(plugInLMM(YsimS[, i], 
                                            predictorLMMmis$fixed.part, predictorLMMmis$random.part, 
                                            predictorLMMmis$reg, predictorLMMmis$con, predictorLMMmis$weights, 
                                           predictorLMMmis$backTrans, predictorLMMmis$thetaFun)$thetaP)
    return(thetaPpluginMis)
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
thetaSim = thetaSim, Ysim = Ysim, errorLMM = errorLMM, errorLMMmis = errorLMMmis))
}
