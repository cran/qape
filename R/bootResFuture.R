bootResFuture <-
function (predictor, B, p, correction) 
{
  N <- nrow(predictor$X)
  YsimF <- function(...) {
    tablsrswrRe <- srswrRe(ranef(predictor$mEst), predictor$reg)$tablsrswrRe
    dfeS <- data.frame(nrw = 1:nrow(predictor$regS), predictor$regS, 
                       eS = predictor$eS)
    dfesamp <- dfeS[sample(nrow(dfeS), N, replace = T), ]
    vsample <- tablsrswrRe$ranef
    names(vsample) <- tablsrswrRe$refNames
    return(predictor$Xbeta + predictor$Z %*% as.matrix(vsample[colnames(predictor$Z)]) + 
             dfesamp$eS)
  }
  YsimT <- function(...) {
    lranefKorekta <- corrRanef(predictor$mEst)
    tablwzlKorekta <- srswrRe(lranefKorekta, predictor$reg)$tablsrswrRe
    eSKorekta <- corrRancomp(predictor$mEst)
    dfeSKorekta <- data.frame(nrw = 1:nrow(predictor$regS), 
                              predictor$regS, eSKorekta)
    dfesampKorekta <- dfeSKorekta[sample(dfeSKorekta$nrw, 
                                         N, replace = T), ]
    vsampleKorekta <- tablwzlKorekta$ranef
    names(vsampleKorekta) <- tablwzlKorekta$refNames
    as.matrix(vsampleKorekta[colnames(predictor$Z)])
    return(predictor$Xbeta + predictor$Z %*% as.matrix(vsampleKorekta[colnames(predictor$Z)]) + 
             dfesampKorekta$eSKorekta)
  }
  if (inherits(predictor, 'EBLUP') == F & inherits(predictor, 'plugInLMM') == F & 
    inherits(predictor, 'ebpLMMne') == F) {
    stop("wrong predictor")
  }
  if (B < 1) {
    stop("B1 must be > 0")
  }
  if (correction == F) {
    Ysim <- matrix(replicate(B, YsimF(predictor)), ncol = B)
  }
  else {
    Ysim <- matrix(replicate(B, YsimT(predictor)), ncol = B)
  }
  YsimS <- matrix(Ysim[predictor$con == 1, ], ncol = B)
  class <- class(predictor)
  
  if (class == "EBLUP") {
    thetaSim <- future_sapply(1:B,  future.seed=TRUE, 
                              future.globals = list(reg=predictor$reg,predictor=predictor, Ysim = Ysim), 
                              function(i) as.numeric(as.vector(predictor$gamma) %*% 
                                                       Ysim[, i]))
    predictorSim <- future_sapply(1:B, future.seed=TRUE, 
                                  future.globals = list(EBLUP=EBLUP, YsimS=YsimS, predictor=predictor, reg=predictor$reg),
                                  function(i) {
                                    thetaPeblup <- as.numeric(EBLUP(YsimS[, i], predictor$fixed.part, 
                                                                    predictor$random.part, predictor$reg, predictor$con, 
                                                                    predictor$gamma, predictor$weights, estMSE = FALSE)$thetaP)
                                    return(thetaPeblup)
                                  })
  }
  

  if (class == "plugInLMM") {
    if (is.null(predictor$backTrans)) {
      predictor$backTrans <- function(x) x
    }
    thetaSim <- future_sapply(1:B, future.seed=TRUE, 
                              future.globals = list(reg=predictor$reg,predictor=predictor, Ysim = Ysim),
                              function(i) as.numeric(predictor$thetaFun(predictor$backTrans(Ysim[, 
                                                                                                 i]))))
    predictorSim <- future_sapply(1:B, future.seed=TRUE, 
                                  future.globals = list(plugInLMM=plugInLMM, YsimS=YsimS, predictor=predictor, reg=predictor$reg),
                                  function(i) {
                                    thetaPplugin <- as.numeric(plugInLMM(YsimS[, i], 
                                                                         predictor$fixed.part, predictor$random.part, 
                                                                         predictor$reg, predictor$con, predictor$weights, 
                                                                         predictor$backTrans, predictor$thetaFun)$thetaP)
                                    return(thetaPplugin)
                                  })
  }
  
if (class == "ebpLMMne") {
    if (is.null(predictor$backTrans)) {
      predictor$backTrans <- function(x) x
    }
    thetaSim <- future_sapply(1:B, future.seed=TRUE, 
                              future.globals = list(reg=predictor$reg,predictor=predictor, Ysim = Ysim),
                              function(i) as.numeric(predictor$thetaFun(predictor$backTrans(Ysim[, 
                                                                                                 i]))))
    predictorSim <- future_sapply(1:B, future.seed=TRUE, 
                                  future.globals = list(ebpLMMne=ebpLMMne, YsimS=YsimS, predictor=predictor, reg=predictor$reg),
                                  function(i) {
                                    thetaPebp <- as.numeric(ebpLMMne(YsimS[, i], predictor$fixed.part, 
                                                                     predictor$division, predictor$reg, predictor$con, 
                                                                     predictor$backTrans, predictor$thetaFun, predictor$L)$thetaP)
                                    return(thetaPebp)
                                  })
}
  
  quantileNaN <- function (x, probs) {
    if (sum(is.nan(x)) > 0) rep(NaN,length(probs)) else {quantile(x, probs)}}
  
  
  error <- matrix((predictorSim - thetaSim), ncol = B)
  return(list(estQAPE = sapply(1:nrow(error), function(i) quantileNaN(abs(error[i, 
  ]), probs = p)), estRMSE = sapply(1:nrow(error), function(i) sqrt((sum(error[i, 
  ]^2))/length(error[i, ]))), predictorSim = predictorSim, 
  thetaSim = thetaSim, Ysim = Ysim, error = error))
}
