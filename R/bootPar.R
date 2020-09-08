bootPar <-
function(predictor, B, p){
  
 if (class(predictor) != 'EBLUP' & class(predictor) != 'plugInLMM'& class(predictor) != 'ebpLMMne') {
      stop("wrong predictor")
    }
 if (B < 1) { 
    stop("B1 must be > 0")
  }
  N <- nrow(predictor$reg)
  
  R <- diag(predictor$sigma2R/predictor$weights, nrow = N, ncol = N)
  
  namesreg <- gsub('\\.', '_', names(predictor$reg))
    
  listVarCorr <- VarCorr(predictor$mEst)
  namesvar <- sub("\\..*", "", names(listVarCorr))
  bmlist <- list()

  for (i in 1:length(listVarCorr)){
      k <- length(unique(predictor$reg[ ,namesvar[i]])) #all v
      bmlist[[i]] <- bdiag(replicate(k,  
                                     as.matrix(as.data.frame(VarCorr(predictor$mEst)[[i]])), simplify = F))
  }  

Gall <- bdiag(bmlist)

Ysim <-  matrix(replicate(B,
              predictor$Xbeta +predictor$Z %*% t(as.matrix(rmvnorm(1, mean = rep(0, ncol(predictor$Z)), 
                                                                    sigma = as.matrix(Gall)))) + 
              t(as.matrix(rmvnorm(1, mean = rep(0, N), sigma = as.matrix(R))))), ncol=B)

YsimS <- Ysim [predictor$con == 1,] 

class <- class(predictor)
if (class == 'EBLUP') {
  thetaSim <- sapply(1:B, function(i) as.numeric(as.vector(predictor$gamma) %*% Ysim[ ,i]))
  predictorSim <- sapply(1:B, function(i) { 
  thetaPeblup <- as.numeric(EBLUP(YsimS[ ,i], predictor$fixed.part, predictor$random.part, predictor$reg,
                                  predictor$con, predictor$gamma, predictor$weights)$thetaP)
  return(thetaPeblup)
  })
}

if (class == 'plugInLMM') {
  if(is.null(predictor$backTrans)){
    predictor$backTrans <- function(x) x
  }
    
  thetaSim <- sapply(1:B, function(i) as.numeric(predictor$thetaFun(predictor$backTrans(Ysim[,i])))) #
 
  predictorSim <- sapply(1:B, function(i) {
    thetaPplugin <- as.numeric(plugInLMM(YsimS[, i], predictor$fixed.part, predictor$random.part, predictor$reg, predictor$con, 
                                       predictor$weights, predictor$backTrans, predictor$thetaFun)$thetaP )
    return(thetaPplugin)
  })
}

if (class == 'ebpLMMne') {
  if(is.null(predictor$backTrans)){
    predictor$backTrans <- function(x) x
  }
  thetaSim <- sapply(1:B, function(i) as.numeric(predictor$thetaFun(predictor$backTrans(Ysim[,i])))) #
  predictorSim <- sapply(1:B, function(i) {
  thetaPebp <- as.numeric(ebpLMMne(YsimS[, i], predictor$fixed.part, predictor$division, predictor$reg, 
                                 predictor$con, predictor$backTrans, predictor$thetaFun, predictor$L)$thetaP)
  return(thetaPebp)
  })
}

error <- matrix((predictorSim - thetaSim),ncol=B) 

return(list(
    estQAPE = sapply(1:nrow(error), function(i) quantile(abs(error[i,]), probs = p)), 
    estRMSE = sapply(1:nrow(error), function(i) sqrt((sum(error[i,]^2))/length(error[i,]))), 
    predictorSim = predictorSim, thetaSim = thetaSim, Ysim = Ysim, error=error 
    ))
}
