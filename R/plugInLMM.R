plugInLMM <-
function(YS, fixed.part, random.part, reg, con, weights, backTrans, thetaFun){
  
  model <- formula(paste('YS', '~', fixed.part, '+', random.part))
  modellm <- formula(paste('YS', '~', fixed.part))
  regS <- subset(reg, con == 1)
  regR <- subset(reg, con == 0)
  
  if(is.null(weights))
      weightS <- rep(1, nrow(regS)) else
      weightS <- subset(weights, con == 1)
  
  Z <- Zfun(model, reg)$Z
  vNames <- Zfun(model, reg)$vNames
  X <- model.matrix(formula(paste('~', fixed.part)), reg)
  
  ZS <- Zfun(model, regS)$Z
  vSNames <- Zfun(model, regS)$vNames
  
  ###estimation
  mEst <- lmer(model, weights = weightS, data.frame(YS, regS, weightS))
  beta <- mEst@beta
  Xbeta <- X %*% beta
  XS <- getME(mEst, name = 'X')
  vS <- as.vector(getME(mEst, name="b"))
  vSDF <- data.frame(vSNames, vS)
  eS <- residuals(mEst)
  
  
  ZR=subset(Z[,(vNames %in% vSNames)], con == 0) 
  
  
  XR <- model.matrix(formula(paste('~', fixed.part)), regR)
  sigma2R <- sigma(mEst)^2
  R <- diag(sigma2R/weightS, nrow = nrow(regS), ncol = nrow(regS))
  G <- sigma(mEst)^2 * crossprod(getME(mEst,"Lambdat")) 
  
  lackOfRanef=(rowSums(Z[,(vNames %in% vSNames)])==0)
  Y=rep(NA,nrow(reg))
  Y[con==1]=YS
  
  if (sum(con==0 &  lackOfRanef==0)>0) {
  Y[con==0 &  lackOfRanef==0]<- predict(mEst, subset(reg, (con == 0 & lackOfRanef ==0)), type = 'response') 
  }
  if (sum(con==0 &  lackOfRanef==1)>0) {
  Y[con==0 &  lackOfRanef==1]<- predict(lm(modellm, weights = weightS, 
                                           data.frame(YS, regS, weightS)), subset(reg, (con == 0 & lackOfRanef ==1)), type = 'response')
  }
  
  YP <- Y[con==0]
  
  if(is.null(backTrans))
    YbackTrans <- Y else
  YbackTrans <- backTrans(Y)

  if(is.null(backTrans))
    YPbackTrans <- YP else
      YPbackTrans <- backTrans(YP)

  thetaP <- thetaFun(YbackTrans)

  outl <- list(thetaP = thetaP, YP = YP, YbackTrans = YbackTrans, YPbackTrans = YPbackTrans,
               beta = beta, Xbeta = Xbeta, sigma2R = sigma2R, R = R, G = G, model = model, 
               mEst = mEst,
               YS = YS, reg = reg, con = con, regS = regS, regR = regR, weights = weightS,
               Z = Z, X = X, ZS = ZS, XR = XR, ZR = ZR, eS = eS, 
               vS = vSDF, fixed.part=fixed.part, 
               random.part=random.part,
               thetaFun = thetaFun, backTrans = backTrans)  
  
  class(outl) = "plugInLMM"
  return(outl)
}
