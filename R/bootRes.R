bootRes <-
function(predictor, B, p, correction) {
    
    if (class(predictor) != 'EBLUP' & class(predictor) != 'plugInLMM'& class(predictor) != 'ebpLMMne') {
      stop("wrong predictor")
    }
    if (B < 1) { 
    stop("B1 must be > 0")
  	}
    N <- nrow(predictor$reg)
    if (correction == F){
      BB = 1
      Ysim <- matrix(data = NA, ncol = B, nrow = N) 

      while(BB <= B){
      tablwzl <- lwzl(ranef(predictor$mEst), predictor$reg)$tablwzl
      
      dfeS <- data.frame(nrw = 1:nrow(predictor$regS),
                   predictor$regS,
                   eS = predictor$eS)
      
      dfesamp=dfeS[sample(nrow(dfeS),N,replace=T),] 
      
      
    
      Ysim[ ,BB] <- predictor$Xbeta + predictor$Z%*%tablwzl$ranef + dfesamp$eS
      BB <- BB + 1
        }
      } 
    else {
      BB = 1
      Ysim <- matrix(data = NA, ncol = B, nrow = N) 

      while(BB <= B){
      lranefKorekta <- corrRanef(predictor$mEst)
      
      tablwzlKorekta <- lwzl(lranefKorekta, predictor$reg)$tablwzl
      eSKorekta <- corrRancomp(predictor$mEst)
      dfeSKorekta <- data.frame(nrw = 1:nrow(predictor$regS), predictor$regS, eSKorekta)
      
      dfesampKorekta=dfeSKorekta[sample(dfeSKorekta$nrw,N,replace=T),] 
      
     Ysim[ ,BB] <- predictor$Xbeta + predictor$Z%*%tablwzlKorekta$ranef + dfesampKorekta$eSKorekta
     BB <- BB + 1
      }
  }
    
    YsimS <- Ysim [predictor$con == 1, ]
    
    class <- class(predictor)
    if (class == 'EBLUP') {
      thetaSim <- sapply(1:B, function(i)
        as.numeric(as.vector(predictor$gamma) %*% Ysim[, i]))
        
      predictorSim <- sapply(1:B, function(i) {
        thetaPeBLUP <- as.numeric(EBLUP(YsimS[ ,i], predictor$fixed.part, predictor$random.part, predictor$reg, predictor$con, predictor$gamma, predictor$weights)$thetaP)
        return(thetaPeBLUP)
      })
    }
    
    if (class == 'plugInLMM') {
      if(is.null(predictor$backTrans)){
        predictor$backTrans <- function(x) x
      }
      
      thetaSim <- sapply(1:B, function(i) as.numeric(predictor$thetaFun(predictor$backTrans(Ysim[,i])))) 
      
      predictorSim <- sapply(1:B, function(i) {
          thetaPplugin<-as.numeric(plugInLMM(YsimS[, i], predictor$fixed.part, predictor$random.part, 
          predictor$reg, predictor$con, predictor$weights, predictor$backTrans, predictor$thetaFun)$thetaP)
        
        return(thetaPplugin)
          
      }) 
    }
    
    if (class == 'ebpLMMne') {
      if(is.null(predictor$backTrans)){
        predictor$backTrans <- function(x) x
      }
      thetaSim <- sapply(1:B, function(i) as.numeric(predictor$thetaFun(predictor$backTrans(Ysim[,i])))) 
      predictorSim <- sapply(1:B, function(i) {
        thetaPebp<-as.numeric(ebpLMMne(YsimS[, i], predictor$fixed.part, predictor$division, predictor$reg, 
                                       predictor$con, predictor$backTrans, predictor$thetaFun, predictor$L)$thetaP )
        return(thetaPebp)
      })
    }
         
    error <- matrix((predictorSim - thetaSim),ncol=B) 
    
    return(list(
      estQAPE = sapply(1:nrow(error), function(i) quantile(abs(error[i,]), probs = p)),
      estRMSE = sapply(1:nrow(error), function(i) sqrt((sum(error[i,]^2))/length(error[i,]))), 
      predictorSim = predictorSim, thetaSim = thetaSim, Ysim = Ysim 
      
    ))
  }
