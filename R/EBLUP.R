EBLUP <-
function(YS,
           fixed.part,
           random.part,
           reg,
           con,
           gamma,
           weights) {
    model <- formula(paste('YS', '~', fixed.part, '+', random.part))
    regS <- subset(reg, con == 1)
    regR <- subset(reg, con == 0)
    gammaS <- subset(gamma, con == 1)
    gammaR <- subset(gamma, con == 0)
    weightS <- subset(weights, con == 1)
    
    Z <- Zfun(model, reg)$Z
    vNames <- Zfun(model, reg)$vNames
    X <- model.matrix(formula(paste('~', fixed.part)), reg)
    
    ZS <- Zfun(model, regS)$Z
    vSNames <- Zfun(model, regS)$vNames
    ZR=subset(Z[,(vNames %in% vSNames)], con == 0) 

    
    mEst <-
      lmer(model, weights = weightS, data.frame(YS, regS, weightS))
    beta <- mEst@beta
    Xbeta <- X %*% beta
    XS <- getME(mEst, name = 'X')
    vS <- as.vector(getME(mEst, name = "b"))
    vSDF <- data.frame(vSNames, vS)
    eS <- residuals(mEst)
    

    XR <- model.matrix(formula(paste('~', fixed.part)), regR)

    thetaP <- as.numeric(gammaS %*% YS + gammaR %*% XR %*% beta +
                              gammaR %*% ZR %*% as.vector(vS))
    sigma2R <- sigma(mEst) ^ 2
    R <- diag(sigma2R / weightS,
              nrow = nrow(regS),
              ncol = nrow(regS))
   G <- sigma(mEst) ^ 2 * crossprod(as.matrix(getME(mEst, "Lambdat")))
    
    outl <-
      list(
        fixed.part=fixed.part,
        random.part=random.part, 
        thetaP = thetaP,
        beta = beta,
        Xbeta = Xbeta,
        sigma2R = sigma2R,
        R = R,
        G = G,
        model = model,
        mEst = mEst,
        YS = YS,
        reg = reg,
        con = con,
        regS = regS,
        regR = regR,
        gamma = gamma,
        gammaS = gammaS,
        gammaR = gammaR,
        weights = weights,
        Z = Z,
        X = X,
        ZS = ZS,
        XR = XR,
        ZR = ZR,
        eS = eS,
        vS = vSDF )
  
    
    
    class(outl) = "EBLUP"
    return(outl)
  }
