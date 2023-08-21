bootRes <-
  function(predictor, B, p, correction) {
    N <- nrow(predictor$X)
    YsimF <- function(...) {
      tablsrswrRe <-
        srswrRe(ranef(predictor$mEst), predictor$reg)$tablsrswrRe
      dfeS <-
        data.frame(nrw = 1:nrow(predictor$regS),
                   predictor$regS,
                   eS = predictor$eS)
      dfesamp <- dfeS[sample(nrow(dfeS), N, replace = T), ]
      
      vsample <- tablsrswrRe$ranef
      names(vsample) <- tablsrswrRe$refNames
      
      #step 1
      return(predictor$Xbeta + predictor$Z %*% as.matrix(vsample[colnames(predictor$Z)]) + dfesamp$eS)
    }
    
    YsimT <- function(...) {
      lranefKorekta <- corrRanef(predictor$mEst)
      tablwzlKorekta <-
        srswrRe(lranefKorekta, predictor$reg)$tablsrswrRe
      eSKorekta <- corrRancomp(predictor$mEst)
      dfeSKorekta <-
        data.frame(nrw = 1:nrow(predictor$regS), predictor$regS, eSKorekta)
      dfesampKorekta <-
        dfeSKorekta[sample(dfeSKorekta$nrw, N, replace = T),]
      
      vsampleKorekta <- tablwzlKorekta$ranef
      names(vsampleKorekta) <- tablwzlKorekta$refNames
      as.matrix(vsampleKorekta[colnames(predictor$Z)])
      
      
      #step 2
      return(
        predictor$Xbeta + predictor$Z %*%  as.matrix(vsampleKorekta[colnames(predictor$Z)]) + dfesampKorekta$eSKorekta
      )
    }
    
    if (inherits(predictor, 'EBLUP') == F &
        inherits(predictor, 'plugInLMM') == F &
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
    
    #step 3
    YsimS <- matrix(Ysim[predictor$con == 1, ], ncol = B)
    
    class <- class(predictor)
    
    #step 4
    if (class == 'EBLUP') {
      thetaSim <- sapply(1:B, function(i)
        as.numeric(as.vector(predictor$gamma) %*% Ysim[, i]))
      
      predictorSim <- sapply(1:B, function(i) {
        thetaPeBLUP <-
          as.numeric(
            EBLUP(
              YsimS[, i],
              predictor$fixed.part,
              predictor$random.part,
              predictor$reg,
              predictor$con,
              predictor$gamma,
              predictor$weights,
              estMSE = FALSE
            )$thetaP
          )
        return(thetaPeBLUP)
      })
    }
    
    if (class == 'plugInLMM') {
      if (is.null(predictor$backTrans)) {
        predictor$backTrans <- function(x)
          x
      }
      
      thetaSim <-
        sapply(1:B, function(i)
          as.numeric(predictor$thetaFun(predictor$backTrans(Ysim[, i])))) #
      
      predictorSim <- sapply(1:B, function(i) {
        thetaPplugin <-
          as.numeric(
            plugInLMM(
              YsimS[, i],
              predictor$fixed.part,
              predictor$random.part,
              predictor$reg,
              predictor$con,
              predictor$weights,
              predictor$backTrans,
              predictor$thetaFun
            )$thetaP
          )
        
        return(thetaPplugin)
        
      })
    }
    
    if (class == 'ebpLMMne') {
      if (is.null(predictor$backTrans)) {
        predictor$backTrans <- function(x)
          x
      }
      thetaSim <-
        sapply(1:B, function(i)
          as.numeric(predictor$thetaFun(predictor$backTrans(Ysim[, i])))) #
      predictorSim <- sapply(1:B, function(i) {
        thetaPebp <-
          as.numeric(
            ebpLMMne(
              YsimS[, i],
              predictor$fixed.part,
              predictor$division,
              predictor$reg,
              predictor$con,
              predictor$backTrans,
              predictor$thetaFun,
              predictor$L
            )$thetaP
          )
        return(thetaPebp)
      })
    }
    
    
    error <- matrix((predictorSim - thetaSim), ncol = B)
    estQAPE <-
      sapply(1:nrow(error), function(i)
        quantileNaN(abs(error[i,]), probs = p))
    estRMSE <-
      sapply(1:nrow(error), function(i)
        sqrt((sum(error[i,] ^ 2)) / length(error[i,])))
    
    error1 <- rbind(estQAPE, estRMSE)  %>% as.data.frame()
    rownames(error1) <- NULL
    error1$names <- c(paste0('QAPE', p), 'estRSME')
    summary <-  melt(error1, id.vars = 'names')
    names(summary) <- c('names', 'thetaFun', '"estAccur')
    
    
    return(
      list(
        estQAPE = estQAPE,
        estRMSE = estRMSE,
        summary = summary,
        predictorSim = predictorSim,
        thetaSim = thetaSim,
        Ysim = Ysim,
        error = error
      )
    )
  }