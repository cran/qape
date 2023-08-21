bootParFutureCor <-
  function(predictor, B, p, ratioR, ratioG)
  {
    if (inherits(predictor, 'EBLUP') == F & inherits(predictor, 'plugInLMM') == F & 
    inherits(predictor, 'ebpLMMne') == F) {
    stop("wrong predictor")
    }
    if (B < 1) {
      stop("B1 must be > 0")
    }
    N <- nrow(predictor$reg)
    R <- diag(predictor$sigma2R / predictor$weights,
              nrow = N,
              ncol = N)
    diag(R) <- diag(R) / ratioR
    namesreg <- gsub("\\.", "_", names(predictor$reg))
    listVarCorr <- VarCorr(predictor$mEst)
    namesvar <- sub("\\..*", "", names(listVarCorr))
    bmlist <- list()
    VarCorrNames <- attributes(listVarCorr)$names
    correct_order <- match(predictor$ZBlockNames, VarCorrNames)
    for (i in 1:length(correct_order)) {
      k <- length(unique(predictor$reg[, namesvar[correct_order[i]]]))
      bmlist[[i]] <-
        bdiag(replicate(k, as.matrix(as.data.frame(
          VarCorr(predictor$mEst)[[correct_order[i]]]
        )),
        simplify = F))
    }
    
    
    Gall <- bdiag(bmlist)
    diag(Gall) <- diag(Gall) / ratioG
    
    if (is.positive.definite(as.matrix(Gall)) == FALSE) {
      positiveDefiniteEstG <- FALSE
      Ysim <-
        as.vector(predictor$Xbeta) + t(chol(as.matrix(R))) %*% matrix(rnorm(N * B), nrow = N)
      cat(
        paste(
          "non-positive definite estimated covariance matrix of random effects - y is generated based on a model without random effects",
          "\n"
        )
      )
    } else
    {
      positiveDefiniteEstG <- TRUE
      Ysim <- as.vector(predictor$Xbeta) +
        predictor$Z %*% (t(chol(as.matrix(Gall))) %*% matrix(rnorm(ncol(predictor$Z) * B), nrow = ncol(predictor$Z))) +
        t(chol(as.matrix(R))) %*% matrix(rnorm(N * B), nrow = N)
    }
    
    
    
    YsimS <- matrix(Ysim[predictor$con == 1,], ncol = B)
    class <- class(predictor)
    if (class == "EBLUP") {
      thetaSim <- future_sapply(1:B,  future.seed = TRUE,
                                future.globals = list(
                                  reg = predictor$reg,
                                  predictor = predictor,
                                  Ysim = Ysim
                                ),
                                function(i)
                                  as.numeric(as.vector(predictor$gamma) %*%
                                               Ysim[, i]))
      predictorSim <- future_sapply(1:B, future.seed = TRUE,
                                    future.globals = list(
                                      EBLUP = EBLUP,
                                      YsimS = YsimS,
                                      predictor = predictor,
                                      reg = predictor$reg
                                    ),
                                    function(i) {
                                      thetaPeblup <- as.numeric(
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
                                      return(thetaPeblup)
                                    })
    }
    if (class == "plugInLMM") {
      if (is.null(predictor$backTrans)) {
        predictor$backTrans <- function(x)
          x
      }
      thetaSim <- future_sapply(1:B, future.seed = TRUE,
                                future.globals = list(
                                  reg = predictor$reg,
                                  predictor = predictor,
                                  Ysim = Ysim
                                ),
                                function(i)
                                  as.numeric(predictor$thetaFun(predictor$backTrans(Ysim[,
                                                                                         i]))))
      predictorSim <- future_sapply(1:B, future.seed = TRUE,
                                    future.globals = list(
                                      plugInLMM = plugInLMM,
                                      YsimS = YsimS,
                                      predictor = predictor,
                                      reg = predictor$reg
                                    ),
                                    function(i) {
                                      thetaPplugin <- as.numeric(
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
    if (class == "ebpLMMne") {
      if (is.null(predictor$backTrans)) {
        predictor$backTrans <- function(x)
          x
      }
      thetaSim <- future_sapply(1:B, future.seed = TRUE,
                                future.globals = list(
                                  reg = predictor$reg,
                                  predictor = predictor,
                                  Ysim = Ysim
                                ),
                                function(i)
                                  as.numeric(predictor$thetaFun(predictor$backTrans(Ysim[,
                                                                                         i]))))
      predictorSim <- future_sapply(1:B, future.seed = TRUE,
                                    future.globals = list(
                                      ebpLMMne = ebpLMMne,
                                      YsimS = YsimS,
                                      predictor = predictor,
                                      reg = predictor$reg
                                    ),
                                    function(i) {
                                      thetaPebp <- as.numeric(
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
      future_sapply(1:nrow(error), function(i)
        quantileNaN(abs(error[i,]), probs = p))
    estRMSE <-
      future_sapply(1:nrow(error), function(i)
        sqrt((sum(error[i,] ^ 2)) / length(error[i,])))
    
    error1 <- rbind(estQAPE, estRMSE) %>% as.data.frame()
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
        error = error,
        positiveDefiniteEstG = positiveDefiniteEstG
      )
    )
  }
