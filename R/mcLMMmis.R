mcLMMmis <-
function(Ypop, predictorLMMmis, predictorLMM, predictorLMM2, K, p, ratioR, ratioG) {
  if (is.null(predictorLMM$backTrans)) {
    predictorLMM$backTrans <- function(x) x
  }
  
  if (is.null(predictorLMM2$backTrans)) {
    predictorLMM2$backTrans <- function(x) x
  }
  
  if (inherits(predictorLMM, 'plugInLMM') == F | inherits(predictorLMMmis, 'plugInLMM') == F | 
    inherits(predictorLMM2, 'plugInLMM') == F) {
    stop("wrong predictor/predictors")
  }
  
  if (!all.equal(predictorLMM$YS, predictorLMM2$YS)){ 
    stop("'YS' object must be the same in both predictors")
  }
  
  if (!all.equal(predictorLMM$reg, predictorLMM2$reg)){ 
    stop("'reg' object must be the same in both predictors")
  }
  
  if (!all.equal(predictorLMM$con,predictorLMM2$con)){ 
    stop("'con' object must be the same in both predictors")
  }
  
  if (!all.equal(predictorLMM$thetaFun, predictorLMM2$thetaFun)){ 
    stop("'thetaFun' object must be the same in both predictors")
  }
  
  if (length(Ypop) != nrow(predictorLMM$reg)){ 
    stop("'Ypop' object must be the population vector")
  }
  
   if (!all.equal((Ypop[predictorLMM$con == 1]), predictorLMM$YS)){ 
    stop("Objects 'Ypop[predictorLMM$con == 1]' and 'predictorLMM$YS' must be the same")
  }
    
  if (!all.equal((Ypop[predictorLMM$con == 1]), predictorLMM2$YS)){ 
    stop("Objects 'Ypop[predictorLMM$con]' and 'predictorLMM2$YS' must be the same")
  }

  
  conPop <- rep(1,length(Ypop))
  
  if (K < 1) {
    stop("K must be > 0")
  }
  

  
  
  mcPop <- bootParFutureCor(plugInLMM(Ypop,predictorLMMmis$fixed.part, 
                                   predictorLMMmis$random.part, predictorLMMmis$reg, 
                                   conPop, predictorLMMmis$weights, 
                                   predictorLMMmis$backTrans, predictorLMMmis$thetaFun), K, p, ratioR, ratioG)
  
  positiveDefiniteEstG <- mcPop$positiveDefiniteEstG
  
  YMc <- mcPop$Ysim
  YMcS <- matrix(YMc[predictorLMMmis$con == 1, ], ncol = K)
  
  thetaMc <- matrix(mcPop$thetaSim, ncol = K)
  EthetaMc <- rowMeans(thetaMc)
  
  
  predictorLMMMc <- future_sapply(1:K, future.seed=TRUE, 
                                  future.globals = list(YMcS = YMcS, plugInLMM = plugInLMM, predictorLMM = predictorLMM, reg = predictorLMM$reg), 
                                  function(i) {
                                    thetaPpluginMc <- as.numeric(plugInLMM(YMcS[, i], 
                                                                           predictorLMM$fixed.part, predictorLMM$random.part,
                                                                           predictorLMM$reg, predictorLMM$con, 
                                                                           predictorLMM$weights, predictorLMM$backTrans,
                                                                           predictorLMM$thetaFun)$thetaP)
                                    return(thetaPpluginMc)
                                  })
  
  predictorLMM2Mc <- future_sapply(1:K, future.seed=TRUE, 
                                 future.globals = list(plugInLMM = plugInLMM, YMcS=YMcS, predictorLMM2 = predictorLMM2, reg = predictorLMM2$reg), 
                                   function(i) {
                                     thetaPpluginLMM2Mc <- as.numeric(plugInLMM(YMcS[, i], 
                                                                                       predictorLMM2$fixed.part, predictorLMM2$random.part,
                                                                                       predictorLMM2$reg, predictorLMM2$con, 
                                                                                       predictorLMM2$weights, predictorLMM2$backTrans,
                                                                                       predictorLMM2$thetaFun)$thetaP)
                                   return(thetaPpluginLMM2Mc)
                                 })
 
  
  errorLMM <- matrix((predictorLMMMc - thetaMc), ncol = K)
  errorLMM2 <- matrix((predictorLMM2Mc - thetaMc), ncol = K)
  QAPElmm <- future_sapply(1:nrow(errorLMM), function(i) quantileNaN(abs(errorLMM[i,]), probs = p)) 
  MSElmm <- future_sapply(1:nrow(errorLMM), function(i) ((sum(errorLMM[i,]^2))/length(errorLMM[i, ])))
  Blmm <- future_sapply(1:nrow(errorLMM), function(i) ((sum(errorLMM[i,]))/length(errorLMM[i, ])))
  rBlmm <- 100 * Blmm / EthetaMc
  RMSElmm <- sqrt(MSElmm)
  rRMSElmm <- 100*RMSElmm/EthetaMc
  QAPElmm2 <- future_sapply(1:nrow(errorLMM2), function(i) quantileNaN(abs(errorLMM2[i,]), probs = p)) 
  MSElmm2 <- future_sapply(1:nrow(errorLMM2), function(i) ((sum(errorLMM2[i,]^2))/length(errorLMM2[i, ])))
  RMSElmm2 <- sqrt(MSElmm2)
  rRMSElmm2 <- 100*RMSElmm2/EthetaMc
  Blmm2 <- future_sapply(1:nrow(errorLMM2), function(i) ((sum(errorLMM2[i,]))/length(errorLMM2[i, ])))
  rBlmm2 <- 100 * Blmm2 / EthetaMc
  
  
  
  return(list(errorLMM = errorLMM,
              errorLMM2 = errorLMM2,
              QAPElmm = QAPElmm,
              RMSElmm = RMSElmm,
              rRMSElmm = rRMSElmm,
              rBlmm = rBlmm,
              QAPElmm2 = QAPElmm2,
              RMSElmm2 = RMSElmm2,
              rRMSElmm2 = rRMSElmm2,
              rBlmm2 = rBlmm2,
              positiveDefiniteEstG = positiveDefiniteEstG))
}
