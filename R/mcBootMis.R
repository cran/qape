mcBootMis <- function(Ypop, predictorLMM, predictorLMMmis, K, B1, B2, p, q) {
  if (is.null(predictorLMM$backTrans)) {
    predictorLMM$backTrans <- function(x) x
  }
  
  if (is.null(predictorLMMmis$backTrans)) {
    predictorLMMmis$backTrans <- function(x) x
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
  
  if (!all.equal(predictorLMM$con,predictorLMMmis$con)){ 
    stop("'con' object must be the same in both predictors")
  }
  
  if (!all.equal(predictorLMM$thetaFun, predictorLMMmis$thetaFun)){ 
    stop("'thetaFun' object must be the same in both predictors")
  }
  
  if (length(Ypop) != nrow(predictorLMM$reg)){ 
    stop("'Ypop' object must be the population vector")
  }
  
  
  if (!all.equal((Ypop[predictorLMM$con == 1]), predictorLMM$backTrans(predictorLMM$YS))){ 
    stop("'Objects 'Ypop[predictorLMM$con == 1]' and 'predictorLMM$backTrans(predictorLMM$YS)' must be the same")
  }
  
  if (!all.equal((Ypop[predictorLMM$con == 1]), predictorLMMmis$backTrans(predictorLMMmis$YS))){ 
    stop("'Objects 'Ypop[predictorLMM$con]' and 'predictorLMMmis$backTrans(predictorLMMmis$YS)' must be the same")
  }
  
  
  conPop <- rep(1,length(Ypop))
  
  if (K < 1) {
    stop("K must be > 0")
  }
  
  if (B1 < 1) {
    stop("B1 must be > 0")
  }
  
  if (B2 < 1) {
    stop("B2 must be > 0")
  }
  

  mcPop <- bootParFuture(plugInLMM(Ypop,predictorLMM$fixed.part, 
                                   predictorLMM$random.part, predictorLMM$reg, 
                                   conPop, predictorLMM$weights, 
                                   predictorLMM$backTrans, predictorLMM$thetaFun), K, p)
  YMc <- mcPop$Ysim
  YMcS <- matrix(YMc[predictorLMM$con == 1, ], ncol = K)
  
 thetaMc <- matrix(mcPop$thetaSim, ncol = K)
 EthetaMc <- rowMeans(thetaMc)
  
  
  predictorLMMMc <- future_sapply(1:K, future.seed=TRUE, 
                                  future.globals = list(YMcS=YMcS, plugInLMM=plugInLMM, predictorLMM=predictorLMM, reg=predictorLMM$reg), 
                                  function(i) {
                                    thetaPpluginMc <- as.numeric(plugInLMM(YMcS[, i], 
                                                                           predictorLMM$fixed.part, predictorLMM$random.part,
                                                                           predictorLMM$reg, predictorLMM$con, 
                                                                           predictorLMM$weights, predictorLMM$backTrans,
                                                                           predictorLMM$thetaFun)$thetaP)
                                    return(thetaPpluginMc)
                                  })
  
  predictorLMMmisMc <- future_sapply(1:K, future.seed=TRUE, 
                                     future.globals = list(plugInLMM=plugInLMM, YMcS=YMcS, predictorLMMmis=predictorLMMmis, reg=predictorLMMmis$reg), 
                                     function(i) {
                                       thetaPpluginLMMmisMc <-  as.numeric(plugInLMM(YMcS[, i], 
                                                                                     predictorLMMmis$fixed.part, predictorLMMmis$random.part,
                                                                                     predictorLMMmis$reg, predictorLMMmis$con, 
                                                                                     predictorLMMmis$weights, predictorLMMmis$backTrans,
                                                                                     predictorLMMmis$thetaFun)$thetaP)
                                       return(thetaPpluginLMMmisMc)
                                     })
  
  quantileNaN <- function (x, probs) {
    if (sum(is.nan(x)) > 0) rep(NaN,length(probs)) else {quantile(x, probs)}}
  
  errorLMM <- matrix((predictorLMMMc - thetaMc), ncol = K)
  errorLMMmis <- matrix((predictorLMMmisMc - thetaMc), ncol = K)
  QAPElmm <- future_sapply(1:nrow(errorLMM), function(i) quantileNaN(abs(errorLMM[i,]), probs = p)) 
  MSElmm <- future_sapply(1:nrow(errorLMM), function(i) ((sum(errorLMM[i,]^2))/length(errorLMM[i, ])))
  Blmm <- future_sapply(1:nrow(errorLMM), function(i) ((sum(errorLMM[i,]))/length(errorLMM[i, ])))
  rBlmm <- 100 * Blmm / EthetaMc
  RMSElmm <- sqrt(MSElmm)
  rRMSElmm <- 100*RMSElmm/EthetaMc
  QAPElmmMis <- future_sapply(1:nrow(errorLMMmis), function(i) quantileNaN(abs(errorLMMmis[i,]), probs = p)) 
  MSElmmMis <- future_sapply(1:nrow(errorLMMmis), function(i) ((sum(errorLMMmis[i,]^2))/length(errorLMMmis[i, ])))
  RMSElmmMis <- sqrt(MSElmmMis)
  rRMSElmmMis <- 100*RMSElmmMis/EthetaMc
  BlmmMis <- future_sapply(1:nrow(errorLMMmis), function(i) ((sum(errorLMMmis[i,]))/length(errorLMMmis[i, ])))
  rBlmmMis <- 100 * BlmmMis / EthetaMc
  
  doubleBootMc <- future_lapply(1:K, future.seed=TRUE, 
                                future.packages = c("matrixcalc", "lme4"), 
                                future.globals = list(B1=B1, B2=B2, p=p, q=q, bootParMis=bootParMis, doubleBootMis=doubleBootMis, YMcS=YMcS, predictorLMM=predictorLMM, reg=predictorLMM$reg, 
                                                      plugInLMM=plugInLMM, predictorLMMmis=predictorLMMmis), function(i) { 
                                                        doubleBootEstAcc <-  doubleBootMis(
                                                          plugInLMM(YMcS[, i], predictorLMM$fixed.part, 
                                                                    predictorLMM$random.part,
                                                                    predictorLMM$reg, predictorLMMmis$con, predictorLMM$weights, 
                                                                    predictorLMM$backTrans, predictorLMM$thetaFun),
                                                          plugInLMM(YMcS[, i], predictorLMMmis$fixed.part, 
                                                                    predictorLMMmis$random.part,
                                                                    predictorLMMmis$reg, predictorLMMmis$con, predictorLMMmis$weights, 
                                                                    predictorLMMmis$backTrans, predictorLMMmis$thetaFun),
                                                          B1,B2, p, q)
                                                        
                                                        return(list(estMSE_param_LMM = doubleBootEstAcc$estMSE_param_LMM,
                                                                    estMSE_db_B2_LMM = doubleBootEstAcc$estMSE_db_B2_LMM,
                                                                    estMSE_db_B2_WDZ_LMM = doubleBootEstAcc$estMSE_db_B2_WDZ_LMM,
                                                                    estMSE_db_B2_HM_LMM = doubleBootEstAcc$estMSE_db_B2_HM_LMM,
                                                                    estMSE_db_1_LMM = doubleBootEstAcc$estMSE_db_1_LMM,
                                                                    estMSE_db_1_WDZ_LMM = doubleBootEstAcc$estMSE_db_1_WDZ_LMM,
                                                                    estMSE_db_1_EF_LMM = doubleBootEstAcc$estMSE_db_1_EF_LMM,
                                                                    estMSE_db_telesc_LMM = doubleBootEstAcc$estMSE_db_telesc_LMM,
                                                                    estMSE_db_telesc_WDZ_LMM = doubleBootEstAcc$estMSE_db_telesc_WDZ_LMM,
                                                                    estMSE_db_telesc_EF_LMM = doubleBootEstAcc$estMSE_db_telesc_EF_LMM,
                                                                    estQAPE_param_LMM = doubleBootEstAcc$estQAPE_param_LMM,
                                                                    estQAPE_db_B2_LMM = doubleBootEstAcc$estQAPE_db_B2_LMM,
                                                                    estQAPE_db_1_LMM = doubleBootEstAcc$estQAPE_db_1_LMM,
                                                                    estQAPE_db_telesc_LMM = doubleBootEstAcc$estQAPE_db_telesc_LMM,
                                                                    estMSE_param_LMMmis = doubleBootEstAcc$estMSE_param_LMMmis,
                                                                    estMSE_db_B2_LMMmis = doubleBootEstAcc$estMSE_db_B2_LMMmis,
                                                                    estMSE_db_B2_WDZ_LMMmis = doubleBootEstAcc$estMSE_db_B2_WDZ_LMMmis,
                                                                    estMSE_db_B2_HM_LMMmis = doubleBootEstAcc$estMSE_db_B2_HM_LMMmis,
                                                                    estMSE_db_1_LMMmis = doubleBootEstAcc$estMSE_db_1_LMMmis,
                                                                    estMSE_db_1_WDZ_LMMmis = doubleBootEstAcc$estMSE_db_1_WDZ_LMMmis,
                                                                    estMSE_db_1_EF_LMMmis = doubleBootEstAcc$estMSE_db_1_EF_LMMmis,
                                                                    estMSE_db_telesc_LMMmis = doubleBootEstAcc$estMSE_db_telesc_LMMmis,
                                                                    estMSE_db_telesc_WDZ_LMMmis = doubleBootEstAcc$estMSE_db_telesc_WDZ_LMMmis,
                                                                    estMSE_db_telesc_EF_LMMmis = doubleBootEstAcc$estMSE_db_telesc_EF_LMMmis,
                                                                    estQAPE_param_LMMmis = doubleBootEstAcc$estQAPE_param_LMMmis,
                                                                    estQAPE_db_B2_LMMmis = doubleBootEstAcc$estQAPE_db_B2_LMMmis,
                                                                    estQAPE_db_1_LMMmis = doubleBootEstAcc$estQAPE_db_1_LMMmis,
                                                                    estQAPE_db_telesc_LMMmis = doubleBootEstAcc$estQAPE_db_telesc_LMMmis,
                                                                    MCpositiveDefiniteEstGlev1 = as.numeric(doubleBootEstAcc$positiveDefiniteEstGlev1),
                                                                    MCpositiveDefiniteEstGlev2 = doubleBootEstAcc$positiveDefiniteEstGlev2))
                                                      })
  
 MCpositiveDefiniteEstGlev1 <- sum(as.numeric(sapply(1:K, function(i) (doubleBootMc[[i]]$MCpositiveDefiniteEstGlev1))))
 MCpositiveDefiniteEstGlev2 <- sum(as.numeric(sapply(1:K, function(i) (doubleBootMc[[i]]$MCpositiveDefiniteEstGlev2))))
  
  
bootResTMc <- future_lapply(1:K, future.seed=TRUE, 
                              future.packages = c("lme4"), 
                              future.globals = list(B1=B1, p=p, bootResMis=bootResMis,YMcS=YMcS, predictorLMM=predictorLMM, reg=predictorLMM$reg, 
                                                    plugInLMM=plugInLMM, predictorLMMmis=predictorLMMmis), function(i) { 
                                                      bootResTEstAcc <-  bootResMis(
                                                        plugInLMM(YMcS[, i], predictorLMM$fixed.part, 
                                                                  predictorLMM$random.part,
                                                                  predictorLMM$reg, predictorLMM$con, predictorLMM$weights, 
                                                                  predictorLMM$backTrans, predictorLMM$thetaFun),
                                                        plugInLMM(YMcS[, i], predictorLMMmis$fixed.part, 
                                                                  predictorLMMmis$random.part,
                                                                  predictorLMMmis$reg, predictorLMMmis$con, predictorLMMmis$weights, 
                                                                  predictorLMMmis$backTrans, predictorLMMmis$thetaFun),
                                                        B1,p, correction = TRUE)
                                                      return(list(
                                                        estQAPElmm = bootResTEstAcc$estQAPElmm, 
                                                        estQAPElmmMis = bootResTEstAcc$estQAPElmmMis,
                                                        estRMSElmm  = bootResTEstAcc$estRMSElmm, 
                                                        estRMSElmmMis = bootResTEstAcc$estRMSElmmMis))    })
  
  
  bootResFMc <- future_lapply(1:K, future.seed=TRUE, 
                              future.packages = c("lme4"), 
                              future.globals = list(B1=B1, p=p, bootResMis=bootResMis,YMcS=YMcS, predictorLMM=predictorLMM, reg=predictorLMM$reg, 
                                                    plugInLMM=plugInLMM, predictorLMMmis=predictorLMMmis), function(i) { 
                                                      bootResFEstAcc <-  bootResMis(
                                                        plugInLMM(YMcS[, i], predictorLMM$fixed.part, 
                                                                  predictorLMM$random.part,
                                                                  predictorLMM$reg, predictorLMM$con, predictorLMM$weights, 
                                                                  predictorLMM$backTrans, predictorLMM$thetaFun),
                                                        plugInLMM(YMcS[, i], predictorLMMmis$fixed.part, 
                                                                  predictorLMMmis$random.part,
                                                                  predictorLMMmis$reg, predictorLMMmis$con, predictorLMMmis$weights, 
                                                                  predictorLMMmis$backTrans, predictorLMMmis$thetaFun),
                                                        B1,p, correction = FALSE)
                                                      return(list(
                                                        estQAPElmm = bootResFEstAcc$estQAPElmm, 
                                                        estQAPElmmMis = bootResFEstAcc$estQAPElmmMis,
                                                        estRMSElmm  = bootResFEstAcc$estRMSElmm, 
                                                        estRMSElmmMis = bootResFEstAcc$estRMSElmmMis))    })
  
  
  MC.estQAPE_rbF_LMM <- lapply(1:K, function(j) {bootResFMc[[j]]$estQAPElmm})
  MC.estQAPE_rbF_LMMmis <- lapply(1:K, function(j) {bootResFMc[[j]]$estQAPElmmMis})
  MC.estRMSE_rbF_LMM <- lapply(1:K, function(j) {bootResFMc[[j]]$estRMSElmm})
  MC.estRMSE_rbF_LMMmis <- lapply(1:K, function(j) {bootResFMc[[j]]$estRMSElmmMis})
  MC.estMSE_rbF_LMM <- lapply(1:K, function(j) {(bootResFMc[[j]]$estRMSElmm)^2})
  MC.estMSE_rbF_LMMmis <- lapply(1:K, function(j) {(bootResFMc[[j]]$estRMSElmmMis)^2})
  
  MC.estQAPE_rbT_LMM <- lapply(1:K, function(j) {bootResTMc[[j]]$estQAPElmm})
  MC.estQAPE_rbT_LMMmis <- lapply(1:K, function(j) {bootResTMc[[j]]$estQAPElmmMis})
  MC.estRMSE_rbT_LMM <- lapply(1:K, function(j) {bootResTMc[[j]]$estRMSElmm})
  MC.estRMSE_rbT_LMMmis <- lapply(1:K, function(j) {bootResTMc[[j]]$estRMSElmmMis})
  MC.estMSE_rbT_LMM <- lapply(1:K, function(j) {(bootResTMc[[j]]$estRMSElmm)^2})
  MC.estMSE_rbT_LMMmis <- lapply(1:K, function(j) {(bootResTMc[[j]]$estRMSElmmMis)^2})
  
  
  MC.estMSE_param_LMM<- lapply(1:K, function(j) {doubleBootMc[[j]]$estMSE_param_LMM})
  MC.estMSE_db_B2_LMM<- lapply(1:K, function(j) {doubleBootMc[[j]]$estMSE_db_B2_LMM})
  MC.estMSE_db_B2_WDZ_LMM<- lapply(1:K, function(j) {doubleBootMc[[j]]$estMSE_db_B2_WDZ_LMM})
  MC.estMSE_db_B2_HM_LMM<- lapply(1:K, function(j) {doubleBootMc[[j]]$estMSE_db_B2_HM_LMM})
  MC.estMSE_db_1_LMM<- lapply(1:K, function(j) {doubleBootMc[[j]]$estMSE_db_1_LMM})
  MC.estMSE_db_1_WDZ_LMM<- lapply(1:K, function(j) {doubleBootMc[[j]]$estMSE_db_1_WDZ_LMM})
  MC.estMSE_db_1_EF_LMM<- lapply(1:K, function(j) {doubleBootMc[[j]]$estMSE_db_1_EF_LMM})
  MC.estMSE_db_telesc_LMM<- lapply(1:K, function(j) {doubleBootMc[[j]]$estMSE_db_telesc_LMM})
  MC.estMSE_db_telesc_WDZ_LMM<- lapply(1:K, function(j) {doubleBootMc[[j]]$estMSE_db_telesc_WDZ_LMM})
  MC.estMSE_db_telesc_EF_LMM<- lapply(1:K, function(j) {doubleBootMc[[j]]$estMSE_db_telesc_EF_LMM})
  
  neg_estMSE_LMM=matrix(c(
    colSums(matrix(as.numeric(unlist(MC.estMSE_param_LMM)<0), byrow=TRUE, nrow = K)),
    colSums(matrix(as.numeric(unlist(MC.estMSE_db_B2_LMM)<0),  byrow=TRUE, nrow = K)),
    colSums(matrix(as.numeric(unlist(MC.estMSE_db_B2_WDZ_LMM)<0),  byrow=TRUE, nrow = K)),
    colSums(matrix(as.numeric(unlist(MC.estMSE_db_B2_HM_LMM)<0),  byrow=TRUE, nrow = K)),
    colSums(matrix(as.numeric(unlist(MC.estMSE_db_1_LMM)<0),  byrow=TRUE, nrow = K)),
    colSums(matrix(as.numeric(unlist(MC.estMSE_db_1_WDZ_LMM)<0),  byrow=TRUE, nrow = K)),
    colSums(matrix(as.numeric(unlist(MC.estMSE_db_1_EF_LMM)<0),  byrow=TRUE, nrow = K)),
    colSums(matrix(as.numeric(unlist(MC.estMSE_db_telesc_LMM)<0),  byrow=TRUE, nrow = K)),
    colSums(matrix(as.numeric(unlist(MC.estMSE_db_telesc_WDZ_LMM)<0),  byrow=TRUE, nrow = K)),
    colSums(matrix(as.numeric(unlist(MC.estMSE_db_telesc_EF_LMM)<0),  byrow=TRUE, nrow = K))), byrow=TRUE, nrow=10)
  names_LMM=(c("estMSE_param_LMM",
               "estMSE_db_B2_LMM",
               "estMSE_db_B2_WDZ_LMM",
               "estMSE_db_B2_HM_LMM",
               "estMSE_db_1_LMM",
               "estMSE_db_1_WDZ_LMM",
               "estMSE_db_1_EF_LMM",
               "estMSE_db_telesc_LMM",
               "estMSE_db_telesc_WDZ_LMM",
               "estMSE_db_telesc_EF_LMM"))
  rownames(neg_estMSE_LMM) <- names_LMM
  
  sqrtNaN <- function (x) {
    a <- rep(NaN,length(x))
    for (i in 1:length(x)){
      if(x[i]>=0) a[i] <- sqrt(x[i])
    }
    c(a)
  }

  MC.estRMSE_param_LMM <- lapply(1:K, function(j) {sqrtNaN(doubleBootMc[[j]]$estMSE_param_LMM)})
  MC.estRMSE_db_B2_LMM <- lapply(1:K, function(j) {sqrtNaN(doubleBootMc[[j]]$estMSE_db_B2_LMM)}) #
  MC.estRMSE_db_B2_WDZ_LMM<- lapply(1:K, function(j) {sqrtNaN(doubleBootMc[[j]]$estMSE_db_B2_WDZ_LMM)}) 
  MC.estRMSE_db_B2_HM_LMM<- lapply(1:K, function(j) {sqrtNaN(doubleBootMc[[j]]$estMSE_db_B2_HM_LMM)})
  MC.estRMSE_db_1_LMM<- lapply(1:K, function(j) {sqrtNaN(doubleBootMc[[j]]$estMSE_db_1_LMM)}) #
  MC.estRMSE_db_1_WDZ_LMM<- lapply(1:K, function(j) {sqrtNaN(doubleBootMc[[j]]$estMSE_db_1_WDZ_LMM)})
  MC.estRMSE_db_1_EF_LMM<- lapply(1:K, function(j) {sqrtNaN(doubleBootMc[[j]]$estMSE_db_1_EF_LMM)}) #
  MC.estRMSE_db_telesc_LMM<- lapply(1:K, function(j) {sqrtNaN(doubleBootMc[[j]]$estMSE_db_telesc_LMM)}) #
  MC.estRMSE_db_telesc_WDZ_LMM<- lapply(1:K, function(j) {sqrtNaN(doubleBootMc[[j]]$estMSE_db_telesc_WDZ_LMM)})
  MC.estRMSE_db_telesc_EF_LMM<- lapply(1:K, function(j) {sqrtNaN(doubleBootMc[[j]]$estMSE_db_telesc_EF_LMM)}) #
  
  MC.estQAPE_param_LMM<- lapply(1:K, function(j) {doubleBootMc[[j]]$estQAPE_param_LMM})
  MC.estQAPE_db_B2_LMM<- lapply(1:K, function(j) {doubleBootMc[[j]]$estQAPE_db_B2_LMM})
  MC.estQAPE_db_1_LMM <- lapply(1:K, function(j) {doubleBootMc[[j]]$estQAPE_db_1_LMM})
  MC.estQAPE_db_telesc_LMM<- lapply(1:K, function(j) {doubleBootMc[[j]]$estQAPE_db_telesc_LMM})
  

  MC.estMSE_param_LMMmis<- lapply(1:K, function(j) {doubleBootMc[[j]]$estMSE_param_LMMmis})
  MC.estMSE_db_B2_LMMmis<- lapply(1:K, function(j) {doubleBootMc[[j]]$estMSE_db_B2_LMMmis})
  MC.estMSE_db_B2_WDZ_LMMmis <- lapply(1:K, function(j) {doubleBootMc[[j]]$estMSE_db_B2_WDZ_LMMmis})
  MC.estMSE_db_B2_HM_LMMmis <- lapply(1:K, function(j) {doubleBootMc[[j]]$estMSE_db_B2_HM_LMMmis})
  MC.estMSE_db_1_LMMmis <- lapply(1:K, function(j) {doubleBootMc[[j]]$estMSE_db_1_LMMmis})
  MC.estMSE_db_1_WDZ_LMMmis<- lapply(1:K, function(j) {doubleBootMc[[j]]$estMSE_db_1_WDZ_LMMmis})
  MC.estMSE_db_1_EF_LMMmis<- lapply(1:K, function(j) {doubleBootMc[[j]]$estMSE_db_1_EF_LMMmis})
  MC.estMSE_db_telesc_LMMmis<- lapply(1:K, function(j) {doubleBootMc[[j]]$estMSE_db_telesc_LMMmis})
  MC.estMSE_db_telesc_WDZ_LMMmis<- lapply(1:K, function(j) {doubleBootMc[[j]]$estMSE_db_telesc_WDZ_LMMmis})
  MC.estMSE_db_telesc_EF_LMMmis<- lapply(1:K, function(j) {doubleBootMc[[j]]$estMSE_db_telesc_EF_LMMmis})
  
  neg_estMSE_LMMmis=matrix(c(
    colSums(matrix(as.numeric(unlist(MC.estMSE_param_LMMmis)<0), byrow=TRUE, nrow = K)),
    colSums(matrix(as.numeric(unlist(MC.estMSE_db_B2_LMMmis)<0),  byrow=TRUE, nrow = K)),
    colSums(matrix(as.numeric(unlist(MC.estMSE_db_B2_WDZ_LMMmis)<0),  byrow=TRUE, nrow = K)),
    colSums(matrix(as.numeric(unlist(MC.estMSE_db_B2_HM_LMMmis)<0),  byrow=TRUE, nrow = K)),
    colSums(matrix(as.numeric(unlist(MC.estMSE_db_1_LMMmis)<0),  byrow=TRUE, nrow = K)),
    colSums(matrix(as.numeric(unlist(MC.estMSE_db_1_WDZ_LMMmis)<0),  byrow=TRUE, nrow = K)),
    colSums(matrix(as.numeric(unlist(MC.estMSE_db_1_EF_LMMmis)<0),  byrow=TRUE, nrow = K)),
    colSums(matrix(as.numeric(unlist(MC.estMSE_db_telesc_LMMmis)<0),  byrow=TRUE, nrow = K)),
    colSums(matrix(as.numeric(unlist(MC.estMSE_db_telesc_WDZ_LMMmis)<0),  byrow=TRUE, nrow = K)),
    colSums(matrix(as.numeric(unlist(MC.estMSE_db_telesc_EF_LMMmis)<0),  byrow=TRUE, nrow = K))), byrow=TRUE, nrow=10)
  names_LMMmis=(c("estMSE_param_LMMmis",
                  "estMSE_db_B2_LMMmis",
                  "estMSE_db_B2_WDZ_LMMmis",
                  "estMSE_db_B2_HM_LMMmis",
                  "estMSE_db_1_LMMmis",
                  "estMSE_db_1_WDZ_LMMmis",
                  "estMSE_db_1_EF_LMMmis",
                  "estMSE_db_telesc_LMMmis",
                  "estMSE_db_telesc_WDZ_LMMmis",
                  "estMSE_db_telesc_EF_LMMmis"))
  rownames(neg_estMSE_LMMmis) <- names_LMMmis
  
  MC.estRMSE_param_LMMmis <- lapply(1:K, function(j) {sqrtNaN(doubleBootMc[[j]]$estMSE_param_LMMmis)})
  MC.estRMSE_db_B2_LMMmis <- lapply(1:K, function(j) {sqrtNaN(doubleBootMc[[j]]$estMSE_db_B2_LMMmis)}) #
  MC.estRMSE_db_B2_WDZ_LMMmis<- lapply(1:K, function(j) {sqrtNaN(doubleBootMc[[j]]$estMSE_db_B2_WDZ_LMMmis)}) 
  MC.estRMSE_db_B2_HM_LMMmis<- lapply(1:K, function(j) {sqrtNaN(doubleBootMc[[j]]$estMSE_db_B2_HM_LMMmis)})
  MC.estRMSE_db_1_LMMmis<- lapply(1:K, function(j) {sqrtNaN(doubleBootMc[[j]]$estMSE_db_1_LMMmis)}) #
  MC.estRMSE_db_1_WDZ_LMMmis<- lapply(1:K, function(j) {sqrtNaN(doubleBootMc[[j]]$estMSE_db_1_WDZ_LMMmis)})
  MC.estRMSE_db_1_EF_LMMmis<- lapply(1:K, function(j) {sqrtNaN(doubleBootMc[[j]]$estMSE_db_1_EF_LMMmis)}) #
  MC.estRMSE_db_telesc_LMMmis<- lapply(1:K, function(j) {sqrtNaN(doubleBootMc[[j]]$estMSE_db_telesc_LMMmis)}) #
  MC.estRMSE_db_telesc_WDZ_LMMmis<- lapply(1:K, function(j) {sqrtNaN(doubleBootMc[[j]]$estMSE_db_telesc_WDZ_LMMmis)})
  MC.estRMSE_db_telesc_EF_LMMmis<- lapply(1:K, function(j) {sqrtNaN(doubleBootMc[[j]]$estMSE_db_telesc_EF_LMMmis)}) #
  
    
  MC.estQAPE_param_LMMmis<- lapply(1:K, function(j) {doubleBootMc[[j]]$estQAPE_param_LMMmis})
  MC.estQAPE_db_B2_LMMmis<- lapply(1:K, function(j) {doubleBootMc[[j]]$estQAPE_db_B2_LMMmis})
  MC.estQAPE_db_1_LMMmis<- lapply(1:K, function(j) {doubleBootMc[[j]]$estQAPE_db_1_LMMmis})
  MC.estQAPE_db_telesc_LMMmis<- lapply(1:K, function(j) {doubleBootMc[[j]]$estQAPE_db_telesc_LMMmis})
 
  
  
  err.estQAPE_rbF_LMM <- lapply(1:K, function(j) {MC.estQAPE_rbF_LMM[[j]] - QAPElmm})
  err.estQAPE_rbF_LMMmis <- lapply(1:K, function(j) {MC.estQAPE_rbF_LMMmis[[j]] - QAPElmmMis})
  err.estRMSE_rbF_LMM <- lapply(1:K, function(j) {MC.estRMSE_rbF_LMM[[j]] - RMSElmm})
  err.estRMSE_rbF_LMMmis <- lapply(1:K, function(j) {MC.estRMSE_rbF_LMMmis[[j]] - RMSElmmMis})
  err.estMSE_rbF_LMM <- lapply(1:K, function(j) {MC.estMSE_rbF_LMM[[j]] - MSElmm})
  err.estMSE_rbF_LMMmis <- lapply(1:K, function(j) {MC.estMSE_rbF_LMMmis[[j]]  - MSElmmMis})
  
  err.estQAPE_rbT_LMM <- lapply(1:K, function(j) {MC.estQAPE_rbT_LMM[[j]]- QAPElmm})
  err.estQAPE_rbT_LMMmis <- lapply(1:K, function(j) {MC.estQAPE_rbT_LMMmis[[j]]- QAPElmmMis})
  err.estRMSE_rbT_LMM <- lapply(1:K, function(j) {MC.estRMSE_rbT_LMM[[j]] - RMSElmm})
  err.estRMSE_rbT_LMMmis <- lapply(1:K, function(j) {MC.estRMSE_rbT_LMMmis[[j]] - RMSElmmMis})
  err.estMSE_rbT_LMM <- lapply(1:K, function(j) {MC.estMSE_rbT_LMM[[j]]  - MSElmm})
  err.estMSE_rbT_LMMmis <- lapply(1:K, function(j) {MC.estMSE_rbT_LMMmis[[j]]  - MSElmmMis})
  
  number_of_predictors <- nrow(thetaMc)
  
  
  rB.estRMSE_rbF_LMM <- 100*colMeans(matrix(unlist(err.estRMSE_rbF_LMM),ncol=number_of_predictors,byrow=TRUE))/RMSElmm
  rRMSE.estRMSE_rbF_LMM <- 100*sqrt(colMeans(matrix(unlist(err.estRMSE_rbF_LMM)^2,ncol=number_of_predictors,byrow=TRUE)))/RMSElmm
  rB.estRMSE_rbF_LMMmis <- 100*colMeans(matrix(unlist(err.estRMSE_rbF_LMMmis),ncol=number_of_predictors,byrow=TRUE))/RMSElmmMis
  rRMSE.estRMSE_rbF_LMMmis <- 100*sqrt(colMeans(matrix(unlist(err.estRMSE_rbF_LMMmis)^2,ncol=number_of_predictors,byrow=TRUE)))/RMSElmmMis
  
  rB.estMSE_rbF_LMM <- 100*colMeans(matrix(unlist(err.estMSE_rbF_LMM),ncol=number_of_predictors,byrow=TRUE))/MSElmm
  rRMSE.estMSE_rbF_LMM <- 100*sqrt(colMeans(matrix(unlist(err.estMSE_rbF_LMM)^2,ncol=number_of_predictors,byrow=TRUE)))/MSElmm
  rB.estMSE_rbF_LMMmis <- 100*colMeans(matrix(unlist(err.estMSE_rbF_LMMmis),ncol=number_of_predictors,byrow=TRUE))/MSElmmMis
  rRMSE.estMSE_rbF_LMMmis <- 100*sqrt(colMeans(matrix(unlist(err.estMSE_rbF_LMMmis)^2,ncol=number_of_predictors,byrow=TRUE)))/MSElmmMis
  
  rB.estQAPE_rbF_LMM <- sapply(1:number_of_predictors, function(j) {
    sapply(1:length(p), function(i) {
      100*mean(sapply(1:K, function(k) {err.estQAPE_rbF_LMM[[k]][i,j]}))/QAPElmm[i,j]
    })
  })
  rownames(rB.estQAPE_rbF_LMM) <- p
  
  rRMSE.estQAPE_rbF_LMM <- sapply(1:number_of_predictors, function(j) {
    sapply(1:length(p), function(i) {
      100*sqrt(mean(sapply(1:K, function(k) {err.estQAPE_rbF_LMM[[k]][i,j]^2})))/QAPElmm[i,j]
    })
  })
  rownames(rRMSE.estQAPE_rbF_LMM) <- p
  
  rB.estQAPE_rbF_LMMmis <- sapply(1:number_of_predictors, function(j) {
    sapply(1:length(p), function(i) {
      100*mean(sapply(1:K, function(k) {err.estQAPE_rbF_LMMmis[[k]][i,j]}))/QAPElmmMis[i,j]
    })
  })
  rownames(rB.estQAPE_rbF_LMMmis) <- p
  
  rRMSE.estQAPE_rbF_LMMmis <- sapply(1:number_of_predictors, function(j) {
    sapply(1:length(p), function(i) {
      100*sqrt(mean(sapply(1:K, function(k) {err.estQAPE_rbF_LMMmis[[k]][i,j]^2})))/QAPElmmMis[i,j]
    })
  })
  rownames(rRMSE.estQAPE_rbF_LMMmis) <- p
  
  
  
  rB.estRMSE_rbT_LMM <- 100*colMeans(matrix(unlist(err.estRMSE_rbT_LMM),ncol=number_of_predictors,byrow=TRUE))/RMSElmm
  rRMSE.estRMSE_rbT_LMM <- 100*sqrt(colMeans(matrix(unlist(err.estRMSE_rbT_LMM)^2,ncol=number_of_predictors,byrow=TRUE)))/RMSElmm
  rB.estRMSE_rbT_LMMmis <- 100*colMeans(matrix(unlist(err.estRMSE_rbT_LMMmis),ncol=number_of_predictors,byrow=TRUE))/RMSElmmMis
  rRMSE.estRMSE_rbT_LMMmis <- 100*sqrt(colMeans(matrix(unlist(err.estRMSE_rbT_LMMmis)^2,ncol=number_of_predictors,byrow=TRUE)))/RMSElmmMis
  
  rB.estMSE_rbT_LMM <- 100*colMeans(matrix(unlist(err.estMSE_rbT_LMM),ncol=number_of_predictors,byrow=TRUE))/MSElmm
  rRMSE.estMSE_rbT_LMM <- 100*sqrt(colMeans(matrix(unlist(err.estMSE_rbT_LMM)^2,ncol=number_of_predictors,byrow=TRUE)))/MSElmm
  rB.estMSE_rbT_LMMmis <- 100*colMeans(matrix(unlist(err.estMSE_rbT_LMMmis),ncol=number_of_predictors,byrow=TRUE))/MSElmmMis
  rRMSE.estMSE_rbT_LMMmis <- 100*sqrt(colMeans(matrix(unlist(err.estMSE_rbT_LMMmis)^2,ncol=number_of_predictors,byrow=TRUE)))/MSElmmMis
  
  rB.estQAPE_rbT_LMM <- sapply(1:number_of_predictors, function(j) {
    sapply(1:length(p), function(i) {
      100*mean(sapply(1:K, function(k) {err.estQAPE_rbT_LMM[[k]][i,j]}))/QAPElmm[i,j]
    })
  })
  rownames(rB.estQAPE_rbT_LMM) <- p
  
  rRMSE.estQAPE_rbT_LMM <- sapply(1:number_of_predictors, function(j) {
    sapply(1:length(p), function(i) {
      100*sqrt(mean(sapply(1:K, function(k) {err.estQAPE_rbT_LMM[[k]][i,j]^2})))/QAPElmm[i,j]
    })
  })
  rownames(rRMSE.estQAPE_rbT_LMM) <- p
  
  rB.estQAPE_rbT_LMMmis <- sapply(1:number_of_predictors, function(j) {
    sapply(1:length(p), function(i) {
      100*mean(sapply(1:K, function(k) {err.estQAPE_rbT_LMMmis[[k]][i,j]}))/QAPElmmMis[i,j]
    })
  })
  rownames(rB.estQAPE_rbT_LMMmis) <- p
  
  rRMSE.estQAPE_rbT_LMMmis <- sapply(1:number_of_predictors, function(j) {
    sapply(1:length(p), function(i) {
      100*sqrt(mean(sapply(1:K, function(k) {err.estQAPE_rbT_LMMmis[[k]][i,j]^2})))/QAPElmmMis[i,j]
    })
  })
  rownames(rRMSE.estQAPE_rbT_LMMmis) <- p
  
 
  err.estMSE_param_LMM <- lapply(1:K, function(j) {MC.estMSE_param_LMM[[j]] - MSElmm})
  err.estMSE_db_B2_LMM <- lapply(1:K, function(j) {MC.estMSE_db_B2_LMM[[j]] - MSElmm})
  err.estMSE_db_B2_WDZ_LMM <- lapply(1:K, function(j) {MC.estMSE_db_B2_WDZ_LMM[[j]] - MSElmm})
  err.estMSE_db_B2_HM_LMM <- lapply(1:K, function(j) {MC.estMSE_db_B2_HM_LMM[[j]] - MSElmm})
  err.estMSE_db_1_LMM <- lapply(1:K, function(j) {MC.estMSE_db_1_LMM[[j]] - MSElmm})
  err.estMSE_db_1_WDZ_LMM <- lapply(1:K, function(j) {MC.estMSE_db_1_WDZ_LMM[[j]] - MSElmm})
  err.estMSE_db_1_EF_LMM <- lapply(1:K, function(j) {MC.estMSE_db_1_EF_LMM[[j]] - MSElmm})
  err.estMSE_db_telesc_LMM <- lapply(1:K, function(j) {MC.estMSE_db_telesc_LMM[[j]] - MSElmm})
  err.estMSE_db_telesc_WDZ_LMM <- lapply(1:K, function(j) {MC.estMSE_db_telesc_WDZ_LMM[[j]] - MSElmm})
  err.estMSE_db_telesc_EF_LMM <- lapply(1:K, function(j) {MC.estMSE_db_telesc_EF_LMM [[j]] - MSElmm})
  
  rB.estMSE_param_LMM <- 100*colMeans(matrix(unlist(err.estMSE_param_LMM),ncol=number_of_predictors,byrow=TRUE))/MSElmm
  rRMSE.estMSE_param_LMM <- 100*sqrt(colMeans(matrix(unlist(err.estMSE_param_LMM)^2,ncol=number_of_predictors,byrow=TRUE)))/MSElmm
  rB.estMSE_db_B2_LMM <- 100*colMeans(matrix(unlist(err.estMSE_db_B2_LMM),ncol=number_of_predictors,byrow=TRUE))/MSElmm
  rRMSE.estMSE_db_B2_LMM <- 100*sqrt(colMeans(matrix(unlist(err.estMSE_db_B2_LMM)^2,ncol=number_of_predictors,byrow=TRUE)))/MSElmm
  rB.estMSE_db_B2_WDZ_LMM <- 100*colMeans(matrix(unlist(err.estMSE_db_B2_WDZ_LMM),ncol=number_of_predictors,byrow=TRUE))/MSElmm
  rRMSE.estMSE_db_B2_WDZ_LMM <- 100*sqrt(colMeans(matrix(unlist(err.estMSE_db_B2_WDZ_LMM)^2,ncol=number_of_predictors,byrow=TRUE)))/MSElmm
  rB.estMSE_db_B2_HM_LMM <- 100*colMeans(matrix(unlist(err.estMSE_db_B2_HM_LMM),ncol=number_of_predictors,byrow=TRUE))/MSElmm
  rRMSE.estMSE_db_B2_HM_LMM <- 100*sqrt(colMeans(matrix(unlist(err.estMSE_db_B2_HM_LMM)^2,ncol=number_of_predictors,byrow=TRUE)))/MSElmm
  rB.estMSE_db_1_LMM <- 100*colMeans(matrix(unlist(err.estMSE_db_1_LMM),ncol=number_of_predictors,byrow=TRUE))/MSElmm
  rRMSE.estMSE_db_1_LMM <- 100*sqrt(colMeans(matrix(unlist(err.estMSE_db_1_LMM)^2,ncol=number_of_predictors,byrow=TRUE)))/MSElmm
  rB.estMSE_db_1_WDZ_LMM <- 100*colMeans(matrix(unlist(err.estMSE_db_1_WDZ_LMM),ncol=number_of_predictors,byrow=TRUE))/MSElmm
  rRMSE.estMSE_db_1_WDZ_LMM <- 100*sqrt(colMeans(matrix(unlist(err.estMSE_db_1_WDZ_LMM)^2,ncol=number_of_predictors,byrow=TRUE)))/MSElmm
  rB.estMSE_db_1_EF_LMM <- 100*colMeans(matrix(unlist(err.estMSE_db_1_EF_LMM),ncol=number_of_predictors,byrow=TRUE))/MSElmm
  rRMSE.estMSE_db_1_EF_LMM <- 100*sqrt(colMeans(matrix(unlist(err.estMSE_db_1_EF_LMM)^2,ncol=number_of_predictors,byrow=TRUE)))/MSElmm
  rB.estMSE_db_telesc_LMM <- 100*colMeans(matrix(unlist(err.estMSE_db_telesc_LMM),ncol=number_of_predictors,byrow=TRUE))/MSElmm
  rRMSE.estMSE_db_telesc_LMM <- 100*sqrt(colMeans(matrix(unlist(err.estMSE_db_telesc_LMM)^2,ncol=number_of_predictors,byrow=TRUE)))/MSElmm
  rB.estMSE_db_telesc_WDZ_LMM <- 100*colMeans(matrix(unlist(err.estMSE_db_telesc_WDZ_LMM),ncol=number_of_predictors,byrow=TRUE))/MSElmm
  rRMSE.estMSE_db_telesc_WDZ_LMM <- 100*sqrt(colMeans(matrix(unlist(err.estMSE_db_telesc_WDZ_LMM)^2,ncol=number_of_predictors,byrow=TRUE)))/MSElmm
  rB.estMSE_db_telesc_EF_LMM <- 100*colMeans(matrix(unlist(err.estMSE_db_telesc_EF_LMM),ncol=number_of_predictors,byrow=TRUE))/MSElmm
  rRMSE.estMSE_db_telesc_EF_LMM <- 100*sqrt(colMeans(matrix(unlist(err.estMSE_db_telesc_EF_LMM)^2,ncol=number_of_predictors,byrow=TRUE)))/MSElmm
  
  err.estRMSE_param_LMM <- lapply(1:K, function(j) {MC.estRMSE_param_LMM[[j]] - RMSElmm})
  err.estRMSE_db_B2_LMM <- lapply(1:K, function(j) {MC.estRMSE_db_B2_LMM[[j]] - RMSElmm})
  err.estRMSE_db_B2_WDZ_LMM <- lapply(1:K, function(j) {MC.estRMSE_db_B2_WDZ_LMM[[j]] - RMSElmm})
  err.estRMSE_db_B2_HM_LMM <- lapply(1:K, function(j) {MC.estRMSE_db_B2_HM_LMM[[j]] - RMSElmm})
  err.estRMSE_db_1_LMM <- lapply(1:K, function(j) {MC.estRMSE_db_1_LMM[[j]] - RMSElmm})
  err.estRMSE_db_1_WDZ_LMM <- lapply(1:K, function(j) {MC.estRMSE_db_1_WDZ_LMM[[j]] - RMSElmm})
  err.estRMSE_db_1_EF_LMM <- lapply(1:K, function(j) {MC.estRMSE_db_1_EF_LMM[[j]] - RMSElmm})
  err.estRMSE_db_telesc_LMM <- lapply(1:K, function(j) {MC.estRMSE_db_telesc_LMM[[j]] - RMSElmm})
  err.estRMSE_db_telesc_WDZ_LMM <- lapply(1:K, function(j) {MC.estRMSE_db_telesc_WDZ_LMM[[j]] - RMSElmm})
  err.estRMSE_db_telesc_EF_LMM <- lapply(1:K, function(j) {MC.estRMSE_db_telesc_EF_LMM [[j]] - RMSElmm})
  
  rB.estRMSE_param_LMM <- 100*colMeans(matrix(unlist(err.estRMSE_param_LMM),ncol=number_of_predictors,byrow=TRUE))/RMSElmm
  rRMSE.estRMSE_param_LMM <- 100*sqrt(colMeans(matrix(unlist(err.estRMSE_param_LMM)^2,ncol=number_of_predictors,byrow=TRUE)))/RMSElmm
  rB.estRMSE_db_B2_LMM <- 100*colMeans(matrix(unlist(err.estRMSE_db_B2_LMM),ncol=number_of_predictors,byrow=TRUE), na.rm = TRUE)/RMSElmm
  rRMSE.estRMSE_db_B2_LMM <- 100*sqrt(colMeans(matrix(unlist(err.estRMSE_db_B2_LMM)^2,ncol=number_of_predictors,byrow=TRUE),na.rm = TRUE))/RMSElmm
  rB.estRMSE_db_B2_WDZ_LMM <- 100*colMeans(matrix(unlist(err.estRMSE_db_B2_WDZ_LMM),ncol=number_of_predictors,byrow=TRUE), na.rm = TRUE)/RMSElmm
  rRMSE.estRMSE_db_B2_WDZ_LMM <- 100*sqrt(colMeans(matrix(unlist(err.estRMSE_db_B2_WDZ_LMM)^2,ncol=number_of_predictors,byrow=TRUE), na.rm = TRUE))/RMSElmm
  rB.estRMSE_db_B2_HM_LMM <- 100*colMeans(matrix(unlist(err.estRMSE_db_B2_HM_LMM),ncol=number_of_predictors,byrow=TRUE), na.rm = TRUE)/RMSElmm
  rRMSE.estRMSE_db_B2_HM_LMM <- 100*sqrt(colMeans(matrix(unlist(err.estRMSE_db_B2_HM_LMM)^2,ncol=number_of_predictors,byrow=TRUE), na.rm = TRUE))/RMSElmm
  rB.estRMSE_db_1_LMM <- 100*colMeans(matrix(unlist(err.estRMSE_db_1_LMM),ncol=number_of_predictors,byrow=TRUE), na.rm = TRUE)/RMSElmm
  rRMSE.estRMSE_db_1_LMM <- 100*sqrt(colMeans(matrix(unlist(err.estRMSE_db_1_LMM)^2,ncol=number_of_predictors,byrow=TRUE), na.rm = TRUE))/RMSElmm
  rB.estRMSE_db_1_WDZ_LMM <- 100*colMeans(matrix(unlist(err.estRMSE_db_1_WDZ_LMM),ncol=number_of_predictors,byrow=TRUE), na.rm = TRUE)/RMSElmm
  rRMSE.estRMSE_db_1_WDZ_LMM <- 100*sqrt(colMeans(matrix(unlist(err.estRMSE_db_1_WDZ_LMM)^2,ncol=number_of_predictors,byrow=TRUE), na.rm = TRUE))/RMSElmm
  rB.estRMSE_db_1_EF_LMM <- 100*colMeans(matrix(unlist(err.estRMSE_db_1_EF_LMM),ncol=number_of_predictors,byrow=TRUE), na.rm = TRUE)/RMSElmm
  rRMSE.estRMSE_db_1_EF_LMM <- 100*sqrt(colMeans(matrix(unlist(err.estRMSE_db_1_EF_LMM)^2,ncol=number_of_predictors,byrow=TRUE), na.rm = TRUE))/RMSElmm
  rB.estRMSE_db_telesc_LMM <- 100*colMeans(matrix(unlist(err.estRMSE_db_telesc_LMM),ncol=number_of_predictors,byrow=TRUE), na.rm = TRUE)/RMSElmm
  rRMSE.estRMSE_db_telesc_LMM <- 100*sqrt(colMeans(matrix(unlist(err.estRMSE_db_telesc_LMM)^2,ncol=number_of_predictors,byrow=TRUE), na.rm = TRUE))/RMSElmm
  rB.estRMSE_db_telesc_WDZ_LMM <- 100*colMeans(matrix(unlist(err.estRMSE_db_telesc_WDZ_LMM),ncol=number_of_predictors,byrow=TRUE), na.rm = TRUE)/RMSElmm
  rRMSE.estRMSE_db_telesc_WDZ_LMM <- 100*sqrt(colMeans(matrix(unlist(err.estRMSE_db_telesc_WDZ_LMM)^2,ncol=number_of_predictors,byrow=TRUE), na.rm = TRUE))/RMSElmm
  rB.estRMSE_db_telesc_EF_LMM <- 100*colMeans(matrix(unlist(err.estRMSE_db_telesc_EF_LMM),ncol=number_of_predictors,byrow=TRUE), na.rm = TRUE)/RMSElmm
  rRMSE.estRMSE_db_telesc_EF_LMM <- 100*sqrt(colMeans(matrix(unlist(err.estRMSE_db_telesc_EF_LMM)^2,ncol=number_of_predictors,byrow=TRUE), na.rm = TRUE))/RMSElmm
  
  
  err.estQAPE_param_LMM <- lapply(1:K, function(j) {MC.estQAPE_param_LMM[[j]] - QAPElmm})
  
  rB.estQAPE_param_LMM <- sapply(1:number_of_predictors, function(j) {
    sapply(1:length(p), function(i) {
      100*mean(sapply(1:K, function(k) {err.estQAPE_param_LMM[[k]][i,j]}))/QAPElmm[i,j]
    })
  })
  rownames(rB.estQAPE_param_LMM) <- p
  
  rRMSE.estQAPE_param_LMM <- sapply(1:number_of_predictors, function(j) {
    sapply(1:length(p), function(i) {
      100*sqrt(mean(sapply(1:K, function(k) {err.estQAPE_param_LMM[[k]][i,j]^2})))/QAPElmm[i,j]
    })
  })
  rownames(rRMSE.estQAPE_param_LMM) <- p
  
  err.estQAPE_db_B2_LMM <- lapply(1:K, function(j) {MC.estQAPE_db_B2_LMM[[j]] - QAPElmm})
  
  rB.estQAPE_db_B2_LMM  <- sapply(1:number_of_predictors, function(j) {
    sapply(1:length(p), function(i) {
      100*mean(sapply(1:K, function(k) {err.estQAPE_db_B2_LMM[[k]][i,j]}))/QAPElmm[i,j]
    })
  })
  rownames(rB.estQAPE_db_B2_LMM ) <- p
  
  rRMSE.estQAPE_db_B2_LMM  <- sapply(1:number_of_predictors, function(j) {
    sapply(1:length(p), function(i) {
      100*sqrt(mean(sapply(1:K, function(k) {err.estQAPE_db_B2_LMM [[k]][i,j]^2})))/QAPElmm[i,j]
    })
  })
  rownames(rRMSE.estQAPE_db_B2_LMM ) <- p
  
  err.estQAPE_db_1_LMM <- lapply(1:K, function(j) {MC.estQAPE_db_1_LMM[[j]] - QAPElmm})
  
  
  rB.estQAPE_db_1_LMM  <- sapply(1:number_of_predictors, function(j) {
    sapply(1:length(p), function(i) {
      100*mean(sapply(1:K, function(k) {err.estQAPE_db_1_LMM[[k]][i,j]}))/QAPElmm[i,j]
    })
  })
  rownames(rB.estQAPE_db_1_LMM) <- p
  
  rRMSE.estQAPE_db_1_LMM <- sapply(1:number_of_predictors, function(j) {
    sapply(1:length(p), function(i) {
      100*sqrt(mean(sapply(1:K, function(k) {err.estQAPE_db_1_LMM[[k]][i,j]^2})))/QAPElmm[i,j]
    })
  })
  rownames(rRMSE.estQAPE_db_1_LMM) <- p
  
  err.estQAPE_db_telesc_LMM <- lapply(1:K, function(j) {MC.estQAPE_db_telesc_LMM[[j]] - QAPElmm})
  
  rB.estQAPE_db_telesc_LMM  <- sapply(1:number_of_predictors, function(j) {
    sapply(1:length(p), function(i) {
      100*mean(sapply(1:K, function(k) {err.estQAPE_db_telesc_LMM[[k]][i,j]}))/QAPElmm[i,j]
    })
  })
  rownames(rB.estQAPE_db_telesc_LMM) <- p
  
  rRMSE.estQAPE_db_telesc_LMM <- sapply(1:number_of_predictors, function(j) {
    sapply(1:length(p), function(i) {
      100*sqrt(mean(sapply(1:K, function(k) {err.estQAPE_db_telesc_LMM[[k]][i,j]^2})))/QAPElmm[i,j]
    })
  })
  rownames(rRMSE.estQAPE_db_telesc_LMM) <- p
  
  err.estMSE_param_LMMmis <- lapply(1:K, function(j) {MC.estMSE_param_LMMmis[[j]] - MSElmmMis})
  err.estMSE_db_B2_LMMmis <- lapply(1:K, function(j) {MC.estMSE_db_B2_LMMmis[[j]] - MSElmmMis})
  err.estMSE_db_B2_WDZ_LMMmis <- lapply(1:K, function(j) {MC.estMSE_db_B2_WDZ_LMMmis[[j]] - MSElmmMis})
  err.estMSE_db_B2_HM_LMMmis <- lapply(1:K, function(j) {MC.estMSE_db_B2_HM_LMMmis[[j]] - MSElmmMis})
  err.estMSE_db_1_LMMmis <- lapply(1:K, function(j) {MC.estMSE_db_1_LMMmis[[j]] - MSElmmMis})
  err.estMSE_db_1_WDZ_LMMmis <- lapply(1:K, function(j) {MC.estMSE_db_1_WDZ_LMMmis[[j]] - MSElmmMis})
  err.estMSE_db_1_EF_LMMmis <- lapply(1:K, function(j) {MC.estMSE_db_1_EF_LMMmis[[j]] - MSElmmMis})
  err.estMSE_db_telesc_LMMmis <- lapply(1:K, function(j) {MC.estMSE_db_telesc_LMMmis[[j]] - MSElmmMis})
  err.estMSE_db_telesc_WDZ_LMMmis <- lapply(1:K, function(j) {MC.estMSE_db_telesc_WDZ_LMMmis[[j]] - MSElmmMis})
  err.estMSE_db_telesc_EF_LMMmis <- lapply(1:K, function(j) {MC.estMSE_db_telesc_EF_LMMmis [[j]] - MSElmmMis})
  
  rB.estMSE_param_LMMmis <- 100*colMeans(matrix(unlist(err.estMSE_param_LMMmis),ncol=number_of_predictors,byrow=TRUE))/MSElmmMis
  rRMSE.estMSE_param_LMMmis <- 100*sqrt(colMeans(matrix(unlist(err.estMSE_param_LMMmis)^2,ncol=number_of_predictors,byrow=TRUE)))/MSElmmMis
  rB.estMSE_db_B2_LMMmis <- 100*colMeans(matrix(unlist(err.estMSE_db_B2_LMMmis),ncol=number_of_predictors,byrow=TRUE))/MSElmmMis
  rRMSE.estMSE_db_B2_LMMmis <- 100*sqrt(colMeans(matrix(unlist(err.estMSE_db_B2_LMMmis)^2,ncol=number_of_predictors,byrow=TRUE)))/MSElmmMis
  rB.estMSE_db_B2_WDZ_LMMmis <- 100*colMeans(matrix(unlist(err.estMSE_db_B2_WDZ_LMMmis),ncol=number_of_predictors,byrow=TRUE))/MSElmmMis
  rRMSE.estMSE_db_B2_WDZ_LMMmis <- 100*sqrt(colMeans(matrix(unlist(err.estMSE_db_B2_WDZ_LMMmis)^2,ncol=number_of_predictors,byrow=TRUE)))/MSElmmMis
  rB.estMSE_db_B2_HM_LMMmis <- 100*colMeans(matrix(unlist(err.estMSE_db_B2_HM_LMMmis),ncol=number_of_predictors,byrow=TRUE))/MSElmmMis
  rRMSE.estMSE_db_B2_HM_LMMmis <- 100*sqrt(colMeans(matrix(unlist(err.estMSE_db_B2_HM_LMMmis)^2,ncol=number_of_predictors,byrow=TRUE)))/MSElmmMis
  rB.estMSE_db_1_LMMmis <- 100*colMeans(matrix(unlist(err.estMSE_db_1_LMMmis),ncol=number_of_predictors,byrow=TRUE))/MSElmmMis
  rRMSE.estMSE_db_1_LMMmis <- 100*sqrt(colMeans(matrix(unlist(err.estMSE_db_1_LMMmis)^2,ncol=number_of_predictors,byrow=TRUE)))/MSElmmMis
  rB.estMSE_db_1_WDZ_LMMmis <- 100*colMeans(matrix(unlist(err.estMSE_db_1_WDZ_LMMmis),ncol=number_of_predictors,byrow=TRUE))/MSElmmMis
  rRMSE.estMSE_db_1_WDZ_LMMmis <- 100*sqrt(colMeans(matrix(unlist(err.estMSE_db_1_WDZ_LMMmis)^2,ncol=number_of_predictors,byrow=TRUE)))/MSElmmMis
  rB.estMSE_db_1_EF_LMMmis <- 100*colMeans(matrix(unlist(err.estMSE_db_1_EF_LMMmis),ncol=number_of_predictors,byrow=TRUE))/MSElmmMis
  rRMSE.estMSE_db_1_EF_LMMmis <- 100*sqrt(colMeans(matrix(unlist(err.estMSE_db_1_EF_LMMmis)^2,ncol=number_of_predictors,byrow=TRUE)))/MSElmmMis
  rB.estMSE_db_telesc_LMMmis <- 100*colMeans(matrix(unlist(err.estMSE_db_telesc_LMMmis),ncol=number_of_predictors,byrow=TRUE))/MSElmmMis
  rRMSE.estMSE_db_telesc_LMMmis <- 100*sqrt(colMeans(matrix(unlist(err.estMSE_db_telesc_LMMmis)^2,ncol=number_of_predictors,byrow=TRUE)))/MSElmmMis
  rB.estMSE_db_telesc_WDZ_LMMmis <- 100*colMeans(matrix(unlist(err.estMSE_db_telesc_WDZ_LMMmis),ncol=number_of_predictors,byrow=TRUE))/MSElmmMis
  rRMSE.estMSE_db_telesc_WDZ_LMMmis <- 100*sqrt(colMeans(matrix(unlist(err.estMSE_db_telesc_WDZ_LMMmis)^2,ncol=number_of_predictors,byrow=TRUE)))/MSElmmMis
  rB.estMSE_db_telesc_EF_LMMmis <- 100*colMeans(matrix(unlist(err.estMSE_db_telesc_EF_LMMmis),ncol=number_of_predictors,byrow=TRUE))/MSElmmMis
  rRMSE.estMSE_db_telesc_EF_LMMmis <- 100*sqrt(colMeans(matrix(unlist(err.estMSE_db_telesc_EF_LMMmis)^2,ncol=number_of_predictors,byrow=TRUE)))/MSElmmMis
  
  err.estRMSE_param_LMMmis <- lapply(1:K, function(j) {MC.estRMSE_param_LMMmis[[j]] - RMSElmmMis})
  err.estRMSE_db_B2_LMMmis <- lapply(1:K, function(j) {MC.estRMSE_db_B2_LMMmis[[j]] - RMSElmmMis})
  err.estRMSE_db_B2_WDZ_LMMmis <- lapply(1:K, function(j) {MC.estRMSE_db_B2_WDZ_LMMmis[[j]] - RMSElmmMis})
  err.estRMSE_db_B2_HM_LMMmis <- lapply(1:K, function(j) {MC.estRMSE_db_B2_HM_LMMmis[[j]] - RMSElmmMis})
  err.estRMSE_db_1_LMMmis <- lapply(1:K, function(j) {MC.estRMSE_db_1_LMMmis[[j]] - RMSElmmMis})
  err.estRMSE_db_1_WDZ_LMMmis <- lapply(1:K, function(j) {MC.estRMSE_db_1_WDZ_LMMmis[[j]] - RMSElmmMis})
  err.estRMSE_db_1_EF_LMMmis <- lapply(1:K, function(j) {MC.estRMSE_db_1_EF_LMMmis[[j]] - RMSElmmMis})
  err.estRMSE_db_telesc_LMMmis <- lapply(1:K, function(j) {MC.estRMSE_db_telesc_LMMmis[[j]] - RMSElmmMis})
  err.estRMSE_db_telesc_WDZ_LMMmis <- lapply(1:K, function(j) {MC.estRMSE_db_telesc_WDZ_LMMmis[[j]] - RMSElmmMis})
  err.estRMSE_db_telesc_EF_LMMmis <- lapply(1:K, function(j) {MC.estRMSE_db_telesc_EF_LMMmis [[j]] - RMSElmmMis})
  
  rB.estRMSE_param_LMMmis <- 100*colMeans(matrix(unlist(err.estRMSE_param_LMMmis),ncol=number_of_predictors,byrow=TRUE))/RMSElmmMis
  rRMSE.estRMSE_param_LMMmis <- 100*sqrt(colMeans(matrix(unlist(err.estRMSE_param_LMMmis)^2,ncol=number_of_predictors,byrow=TRUE)))/RMSElmmMis
  rB.estRMSE_db_B2_LMMmis <- 100*colMeans(matrix(unlist(err.estRMSE_db_B2_LMMmis),ncol=number_of_predictors,byrow=TRUE), na.rm = TRUE)/RMSElmmMis
  rRMSE.estRMSE_db_B2_LMMmis <- 100*sqrt(colMeans(matrix(unlist(err.estRMSE_db_B2_LMMmis)^2,ncol=number_of_predictors,byrow=TRUE),na.rm = TRUE))/RMSElmmMis
  rB.estRMSE_db_B2_WDZ_LMMmis <- 100*colMeans(matrix(unlist(err.estRMSE_db_B2_WDZ_LMMmis),ncol=number_of_predictors,byrow=TRUE), na.rm = TRUE)/RMSElmmMis
  rRMSE.estRMSE_db_B2_WDZ_LMMmis <- 100*sqrt(colMeans(matrix(unlist(err.estRMSE_db_B2_WDZ_LMMmis)^2,ncol=number_of_predictors,byrow=TRUE), na.rm = TRUE))/RMSElmmMis
  rB.estRMSE_db_B2_HM_LMMmis <- 100*colMeans(matrix(unlist(err.estRMSE_db_B2_HM_LMMmis),ncol=number_of_predictors,byrow=TRUE), na.rm = TRUE)/RMSElmmMis
  rRMSE.estRMSE_db_B2_HM_LMMmis <- 100*sqrt(colMeans(matrix(unlist(err.estRMSE_db_B2_HM_LMMmis)^2,ncol=number_of_predictors,byrow=TRUE), na.rm = TRUE))/RMSElmmMis
  rB.estRMSE_db_1_LMMmis <- 100*colMeans(matrix(unlist(err.estRMSE_db_1_LMMmis),ncol=number_of_predictors,byrow=TRUE), na.rm = TRUE)/RMSElmmMis
  rRMSE.estRMSE_db_1_LMMmis <- 100*sqrt(colMeans(matrix(unlist(err.estRMSE_db_1_LMMmis)^2,ncol=number_of_predictors,byrow=TRUE), na.rm = TRUE))/RMSElmmMis
  rB.estRMSE_db_1_WDZ_LMMmis <- 100*colMeans(matrix(unlist(err.estRMSE_db_1_WDZ_LMMmis),ncol=number_of_predictors,byrow=TRUE), na.rm = TRUE)/RMSElmmMis
  rRMSE.estRMSE_db_1_WDZ_LMMmis <- 100*sqrt(colMeans(matrix(unlist(err.estRMSE_db_1_WDZ_LMMmis)^2,ncol=number_of_predictors,byrow=TRUE), na.rm = TRUE))/RMSElmmMis
  rB.estRMSE_db_1_EF_LMMmis <- 100*colMeans(matrix(unlist(err.estRMSE_db_1_EF_LMMmis),ncol=number_of_predictors,byrow=TRUE), na.rm = TRUE)/RMSElmmMis
  rRMSE.estRMSE_db_1_EF_LMMmis <- 100*sqrt(colMeans(matrix(unlist(err.estRMSE_db_1_EF_LMMmis)^2,ncol=number_of_predictors,byrow=TRUE), na.rm = TRUE))/RMSElmmMis
  rB.estRMSE_db_telesc_LMMmis <- 100*colMeans(matrix(unlist(err.estRMSE_db_telesc_LMMmis),ncol=number_of_predictors,byrow=TRUE), na.rm = TRUE)/RMSElmmMis
  rRMSE.estRMSE_db_telesc_LMMmis <- 100*sqrt(colMeans(matrix(unlist(err.estRMSE_db_telesc_LMMmis)^2,ncol=number_of_predictors,byrow=TRUE), na.rm = TRUE))/RMSElmmMis
  rB.estRMSE_db_telesc_WDZ_LMMmis <- 100*colMeans(matrix(unlist(err.estRMSE_db_telesc_WDZ_LMMmis),ncol=number_of_predictors,byrow=TRUE), na.rm = TRUE)/RMSElmmMis
  rRMSE.estRMSE_db_telesc_WDZ_LMMmis <- 100*sqrt(colMeans(matrix(unlist(err.estRMSE_db_telesc_WDZ_LMMmis)^2,ncol=number_of_predictors,byrow=TRUE), na.rm = TRUE))/RMSElmmMis
  rB.estRMSE_db_telesc_EF_LMMmis <- 100*colMeans(matrix(unlist(err.estRMSE_db_telesc_EF_LMMmis),ncol=number_of_predictors,byrow=TRUE), na.rm = TRUE)/RMSElmmMis
  rRMSE.estRMSE_db_telesc_EF_LMMmis <- 100*sqrt(colMeans(matrix(unlist(err.estRMSE_db_telesc_EF_LMMmis)^2,ncol=number_of_predictors,byrow=TRUE), na.rm = TRUE))/RMSElmmMis
  
  err.estQAPE_param_LMMmis <- lapply(1:K, function(j) {MC.estQAPE_param_LMMmis[[j]] - QAPElmmMis})
  
  rB.estQAPE_param_LMMmis <- sapply(1:number_of_predictors, function(j) {
    sapply(1:length(p), function(i) {
      100*mean(sapply(1:K, function(k) {err.estQAPE_param_LMMmis[[k]][i,j]}))/QAPElmmMis[i,j]
    })
  })
  rownames(rB.estQAPE_param_LMMmis) <- p
  
  rRMSE.estQAPE_param_LMMmis <- sapply(1:number_of_predictors, function(j) {
    sapply(1:length(p), function(i) {
      100*sqrt(mean(sapply(1:K, function(k) {err.estQAPE_param_LMMmis[[k]][i,j]^2})))/QAPElmmMis[i,j]
    })
  })
  rownames(rRMSE.estQAPE_param_LMMmis) <- p
  
  err.estQAPE_db_B2_LMMmis <- lapply(1:K, function(j) {MC.estQAPE_db_B2_LMMmis[[j]] - QAPElmmMis})
  
  rB.estQAPE_db_B2_LMMmis  <- sapply(1:number_of_predictors, function(j) {
    sapply(1:length(p), function(i) {
      100*mean(sapply(1:K, function(k) {err.estQAPE_db_B2_LMMmis[[k]][i,j]}))/QAPElmmMis[i,j]
    })
  })
  rownames(rB.estQAPE_db_B2_LMMmis ) <- p
  
  rRMSE.estQAPE_db_B2_LMMmis  <- sapply(1:number_of_predictors, function(j) {
    sapply(1:length(p), function(i) {
      100*sqrt(mean(sapply(1:K, function(k) {err.estQAPE_db_B2_LMMmis [[k]][i,j]^2})))/QAPElmmMis[i,j]
    })
  })
  rownames(rRMSE.estQAPE_db_B2_LMMmis ) <- p
  
  err.estQAPE_db_1_LMMmis <- lapply(1:K, function(j) {MC.estQAPE_db_1_LMMmis[[j]] - QAPElmmMis})
  
  rB.estQAPE_db_1_LMMmis  <- sapply(1:number_of_predictors, function(j) {
    sapply(1:length(p), function(i) {
      100*mean(sapply(1:K, function(k) {err.estQAPE_db_1_LMMmis[[k]][i,j]}))/QAPElmmMis[i,j]
    })
  })
  rownames(rB.estQAPE_db_1_LMMmis) <- p
  
  rRMSE.estQAPE_db_1_LMMmis <- sapply(1:number_of_predictors, function(j) {
    sapply(1:length(p), function(i) {
      100*sqrt(mean(sapply(1:K, function(k) {err.estQAPE_db_1_LMMmis[[k]][i,j]^2})))/QAPElmmMis[i,j]
    })
  })
  rownames(rRMSE.estQAPE_db_1_LMMmis) <- p
  
  err.estQAPE_db_telesc_LMMmis <- lapply(1:K, function(j) {MC.estQAPE_db_telesc_LMMmis[[j]] - QAPElmmMis})
  
  rB.estQAPE_db_telesc_LMMmis  <- sapply(1:number_of_predictors, function(j) {
    sapply(1:length(p), function(i) {
      100*mean(sapply(1:K, function(k) {err.estQAPE_db_telesc_LMMmis[[k]][i,j]}))/QAPElmmMis[i,j]
    })
  })
  rownames(rB.estQAPE_db_telesc_LMMmis) <- p
  
  rRMSE.estQAPE_db_telesc_LMMmis <- sapply(1:number_of_predictors, function(j) {
    sapply(1:length(p), function(i) {
      100*sqrt(mean(sapply(1:K, function(k) {err.estQAPE_db_telesc_LMMmis[[k]][i,j]^2})))/QAPElmmMis[i,j]
    })
  })
  rownames(rRMSE.estQAPE_db_telesc_LMMmis) <- p
  
  return(list(QAPElmm = QAPElmm,
              RMSElmm = RMSElmm,
              rRMSElmm = rRMSElmm,
              rBlmm = rBlmm,
              QAPElmmMis = QAPElmmMis,
              RMSElmmMis = RMSElmmMis,
              rRMSElmmMis = rRMSElmmMis,
              rBlmmMis = rBlmmMis,
              rB.estRMSE_rbF_LMM = rB.estRMSE_rbF_LMM, 
              rRMSE.estRMSE_rbF_LMM = rRMSE.estRMSE_rbF_LMM,
              rB.estRMSE_rbF_LMMmis = rB.estRMSE_rbF_LMMmis, 
              rRMSE.estRMSE_rbF_LMMmis = rRMSE.estRMSE_rbF_LMMmis, 
              rB.estMSE_rbF_LMM = rB.estMSE_rbF_LMM, 
              rRMSE.estMSE_rbF_LMM = rRMSE.estMSE_rbF_LMM, 
              rB.estMSE_rbF_LMMmis = rB.estMSE_rbF_LMMmis,
              rRMSE.estMSE_rbF_LMMmis = rRMSE.estMSE_rbF_LMMmis, 
              rB.estQAPE_rbF_LMM = rB.estQAPE_rbF_LMM,
              rRMSE.estQAPE_rbF_LMM = rRMSE.estQAPE_rbF_LMM,
              rB.estQAPE_rbF_LMMmis = rB.estQAPE_rbF_LMMmis,
              rRMSE.estQAPE_rbF_LMMmis = rRMSE.estQAPE_rbF_LMMmis,
              rB.estRMSE_rbT_LMM = rB.estRMSE_rbT_LMM, 
              rRMSE.estRMSE_rbT_LMM = rRMSE.estRMSE_rbT_LMM,
              rB.estRMSE_rbT_LMMmis = rB.estRMSE_rbT_LMMmis, 
              rRMSE.estRMSE_rbT_LMMmis = rRMSE.estRMSE_rbT_LMMmis, 
              rB.estMSE_rbT_LMM = rB.estMSE_rbT_LMM, 
              rRMSE.estMSE_rbT_LMM = rRMSE.estMSE_rbT_LMM, 
              rB.estMSE_rbT_LMMmis = rB.estMSE_rbT_LMMmis,
              rRMSE.estMSE_rbT_LMMmis = rRMSE.estMSE_rbT_LMMmis, 
              rB.estQAPE_rbT_LMM = rB.estQAPE_rbT_LMM,
              rRMSE.estQAPE_rbT_LMM = rRMSE.estQAPE_rbT_LMM,
              rB.estQAPE_rbT_LMMmis = rB.estQAPE_rbT_LMMmis,
              rRMSE.estQAPE_rbT_LMMmis = rRMSE.estQAPE_rbT_LMMmis,
              neg_estMSE_LMM = neg_estMSE_LMM, 
              neg_estMSE_LMMmis = neg_estMSE_LMMmis,
              rB.estMSE_param_LMMmis = rB.estMSE_param_LMMmis,
              rRMSE.estMSE_param_LMMmis = rRMSE.estMSE_param_LMMmis,
              rB.estMSE_db_B2_LMMmis = rB.estMSE_db_B2_LMMmis,
              rRMSE.estMSE_db_B2_LMMmis = rRMSE.estMSE_db_B2_LMMmis,
              rB.estMSE_db_B2_WDZ_LMMmis = rB.estMSE_db_B2_WDZ_LMMmis,
              rRMSE.estMSE_db_B2_WDZ_LMMmis = rRMSE.estMSE_db_B2_WDZ_LMMmis,
              rB.estMSE_db_B2_HM_LMMmis = rB.estMSE_db_B2_HM_LMMmis,
              rRMSE.estMSE_db_B2_HM_LMMmis = rRMSE.estMSE_db_B2_HM_LMMmis,
              rB.estMSE_db_1_LMMmis = rB.estMSE_db_1_LMMmis,
              rRMSE.estMSE_db_1_LMMmis = rRMSE.estMSE_db_1_LMMmis,
              rB.estMSE_db_1_WDZ_LMMmis = rB.estMSE_db_1_WDZ_LMMmis,
              rRMSE.estMSE_db_1_WDZ_LMMmis = rRMSE.estMSE_db_1_WDZ_LMMmis,
              rB.estMSE_db_1_EF_LMMmis = rB.estMSE_db_1_EF_LMMmis,
              rRMSE.estMSE_db_1_EF_LMMmis = rRMSE.estMSE_db_1_EF_LMMmis,
              rB.estMSE_db_telesc_LMMmis = rB.estMSE_db_telesc_LMMmis,
              rRMSE.estMSE_db_telesc_LMMmis = rRMSE.estMSE_db_telesc_LMMmis,
              rB.estMSE_db_telesc_WDZ_LMMmis = rB.estMSE_db_telesc_WDZ_LMMmis,
              rRMSE.estMSE_db_telesc_WDZ_LMMmis = rRMSE.estMSE_db_telesc_WDZ_LMMmis,
              rB.estMSE_db_telesc_EF_LMMmis = rB.estMSE_db_telesc_EF_LMMmis,
              rRMSE.estMSE_db_telesc_EF_LMMmis = rRMSE.estMSE_db_telesc_EF_LMMmis,
              rB.estRMSE_param_LMMmis = rB.estRMSE_param_LMMmis,
              rRMSE.estRMSE_param_LMMmis = rRMSE.estRMSE_param_LMMmis,
              rB.estRMSE_db_B2_LMMmis = rB.estRMSE_db_B2_LMMmis,
              rRMSE.estRMSE_db_B2_LMMmis = rRMSE.estRMSE_db_B2_LMMmis,
              rB.estRMSE_db_B2_WDZ_LMMmis = rB.estRMSE_db_B2_WDZ_LMMmis,
              rRMSE.estRMSE_db_B2_WDZ_LMMmis = rRMSE.estRMSE_db_B2_WDZ_LMMmis,
              rB.estRMSE_db_B2_HM_LMMmis = rB.estRMSE_db_B2_HM_LMMmis,
              rRMSE.estRMSE_db_B2_HM_LMMmis = rRMSE.estRMSE_db_B2_HM_LMMmis,
              rB.estRMSE_db_1_LMMmis = rB.estRMSE_db_1_LMMmis,
              rRMSE.estRMSE_db_1_LMMmis = rRMSE.estRMSE_db_1_LMMmis,
              rB.estRMSE_db_1_WDZ_LMMmis = rB.estRMSE_db_1_WDZ_LMMmis,
              rRMSE.estRMSE_db_1_WDZ_LMMmis = rRMSE.estRMSE_db_1_WDZ_LMMmis,
              rB.estRMSE_db_1_EF_LMMmis = rB.estRMSE_db_1_EF_LMMmis,
              rRMSE.estRMSE_db_1_EF_LMMmis = rRMSE.estRMSE_db_1_EF_LMMmis,
              rB.estRMSE_db_telesc_LMMmis = rB.estRMSE_db_telesc_LMMmis,
              rRMSE.estRMSE_db_telesc_LMMmis = rRMSE.estRMSE_db_telesc_LMMmis,
              rB.estRMSE_db_telesc_WDZ_LMMmis = rB.estRMSE_db_telesc_WDZ_LMMmis,
              rRMSE.estRMSE_db_telesc_WDZ_LMMmis = rRMSE.estRMSE_db_telesc_WDZ_LMMmis,
              rB.estRMSE_db_telesc_EF_LMMmis = rB.estRMSE_db_telesc_EF_LMMmis,
              rRMSE.estRMSE_db_telesc_EF_LMMmis = rRMSE.estRMSE_db_telesc_EF_LMMmis,
              rB.estQAPE_param_LMMmis = rB.estQAPE_param_LMMmis,
              rRMSE.estQAPE_param_LMMmis = rRMSE.estQAPE_param_LMMmis,
              rB.estQAPE_db_B2_LMMmis = rB.estQAPE_db_B2_LMMmis,
              rRMSE.estQAPE_db_B2_LMMmis = rRMSE.estQAPE_db_B2_LMMmis,
              rB.estQAPE_db_1_LMMmis = rB.estQAPE_db_1_LMMmis,
              rRMSE.estQAPE_db_1_LMMmis = rRMSE.estQAPE_db_1_LMMmis,
              rB.estQAPE_db_telesc_LMMmis = rB.estQAPE_db_telesc_LMMmis,
              rRMSE.estQAPE_db_telesc_LMMmis = rRMSE.estQAPE_db_telesc_LMMmis,
              rB.estMSE_param_LMM = rB.estMSE_param_LMM,
              rRMSE.estMSE_param_LMM = rRMSE.estMSE_param_LMM,
              rB.estMSE_db_B2_LMM = rB.estMSE_db_B2_LMM,
              rRMSE.estMSE_db_B2_LMM = rRMSE.estMSE_db_B2_LMM,
              rB.estMSE_db_B2_WDZ_LMM = rB.estMSE_db_B2_WDZ_LMM,
              rRMSE.estMSE_db_B2_WDZ_LMM = rRMSE.estMSE_db_B2_WDZ_LMM,
              rB.estMSE_db_B2_HM_LMM = rB.estMSE_db_B2_HM_LMM,
              rRMSE.estMSE_db_B2_HM_LMM = rRMSE.estMSE_db_B2_HM_LMM,
              rB.estMSE_db_1_LMM = rB.estMSE_db_1_LMM,
              rRMSE.estMSE_db_1_LMM = rRMSE.estMSE_db_1_LMM,
              rB.estMSE_db_1_WDZ_LMM = rB.estMSE_db_1_WDZ_LMM,
              rRMSE.estMSE_db_1_WDZ_LMM = rRMSE.estMSE_db_1_WDZ_LMM,
              rB.estMSE_db_1_EF_LMM = rB.estMSE_db_1_EF_LMM,
              rRMSE.estMSE_db_1_EF_LMM = rRMSE.estMSE_db_1_EF_LMM,
              rB.estMSE_db_telesc_LMM = rB.estMSE_db_telesc_LMM,
              rRMSE.estMSE_db_telesc_LMM = rRMSE.estMSE_db_telesc_LMM,
              rB.estMSE_db_telesc_WDZ_LMM = rB.estMSE_db_telesc_WDZ_LMM,
              rRMSE.estMSE_db_telesc_WDZ_LMM = rRMSE.estMSE_db_telesc_WDZ_LMM,
              rB.estMSE_db_telesc_EF_LMM = rB.estMSE_db_telesc_EF_LMM,
              rRMSE.estMSE_db_telesc_EF_LMM = rRMSE.estMSE_db_telesc_EF_LMM,
              rB.estRMSE_param_LMM = rB.estRMSE_param_LMM,
              rRMSE.estRMSE_param_LMM = rRMSE.estRMSE_param_LMM,
              rB.estRMSE_db_B2_LMM = rB.estRMSE_db_B2_LMM,
              rRMSE.estRMSE_db_B2_LMM = rRMSE.estRMSE_db_B2_LMM,
              rB.estRMSE_db_B2_WDZ_LMM = rB.estRMSE_db_B2_WDZ_LMM,
              rRMSE.estRMSE_db_B2_WDZ_LMM = rRMSE.estRMSE_db_B2_WDZ_LMM,
              rB.estRMSE_db_B2_HM_LMM = rB.estRMSE_db_B2_HM_LMM,
              rRMSE.estRMSE_db_B2_HM_LMM = rRMSE.estRMSE_db_B2_HM_LMM,
              rB.estRMSE_db_1_LMM = rB.estRMSE_db_1_LMM,
              rRMSE.estRMSE_db_1_LMM = rRMSE.estRMSE_db_1_LMM,
              rB.estRMSE_db_1_WDZ_LMM = rB.estRMSE_db_1_WDZ_LMM,
              rRMSE.estRMSE_db_1_WDZ_LMM = rRMSE.estRMSE_db_1_WDZ_LMM,
              rB.estRMSE_db_1_EF_LMM = rB.estRMSE_db_1_EF_LMM,
              rRMSE.estRMSE_db_1_EF_LMM = rRMSE.estRMSE_db_1_EF_LMM,
              rB.estRMSE_db_telesc_LMM = rB.estRMSE_db_telesc_LMM,
              rRMSE.estRMSE_db_telesc_LMM = rRMSE.estRMSE_db_telesc_LMM,
              rB.estRMSE_db_telesc_WDZ_LMM = rB.estRMSE_db_telesc_WDZ_LMM,
              rRMSE.estRMSE_db_telesc_WDZ_LMM = rRMSE.estRMSE_db_telesc_WDZ_LMM,
              rB.estRMSE_db_telesc_EF_LMM = rB.estRMSE_db_telesc_EF_LMM,
              rRMSE.estRMSE_db_telesc_EF_LMM = rRMSE.estRMSE_db_telesc_EF_LMM,
              rB.estQAPE_param_LMM = rB.estQAPE_param_LMM,
              rRMSE.estQAPE_param_LMM = rRMSE.estQAPE_param_LMM,
              rB.estQAPE_db_B2_LMM = rB.estQAPE_db_B2_LMM,
              rRMSE.estQAPE_db_B2_LMM = rRMSE.estQAPE_db_B2_LMM,
              rB.estQAPE_db_1_LMM = rB.estQAPE_db_1_LMM,
              rRMSE.estQAPE_db_1_LMM = rRMSE.estQAPE_db_1_LMM,
              rB.estQAPE_db_telesc_LMM = rB.estQAPE_db_telesc_LMM,
              rRMSE.estQAPE_db_telesc_LMM = rRMSE.estQAPE_db_telesc_LMM,
              MCpositiveDefiniteEstGlev1 = MCpositiveDefiniteEstGlev1,
              MCpositiveDefiniteEstGlev2 = MCpositiveDefiniteEstGlev2))
}