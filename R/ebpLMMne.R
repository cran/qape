ebpLMMne <-
  function(YS,
           fixed.part,
           division,
           reg,
           con,
           backTrans,
           thetaFun,
           L)
  {
    N = nrow(reg)
    random.part = paste("(1|", paste(division), ")")
    model <- formula(paste("YS", "~", fixed.part,
                           "+", random.part))
    regS <- subset(reg, con == 1)
    regR <- subset(reg, con == 0)
    weights = rep(1, N)
    mEst <- lmer(model, data.frame(YS, regS))
    Zobj <- Zfun(model, reg)
    Z <- Zobj$Z
    ZBlockNames <- Zobj$ZBlockNames
    vNames <- make.unique(Zfun(model, reg)$vNames,sep = ".") 
    colnames(Z) <- vNames
    X <- model.matrix(formula(paste("~", fixed.part)), reg)
    ZS <- getME(mEst, name = "Z")
    vSNames <- make.unique(colnames(ZS),sep = ".") 
    colnames(ZS) <- vSNames 
    sigmaR = sigma(mEst)
    sigma2R = sigmaR ^ 2
    S2v = as.numeric(VarCorr(mEst)[1])
    Sv = sqrt(S2v)
    beta <- mEst@beta
    Xbeta <- X %*% beta
    XS <- getME(mEst, name = "X")
    vS <- as.vector(getME(mEst, name = "b"))
    vSDF <- data.frame(vSNames, vS)
    eS <- residuals(mEst)
    ZR <- Z[(con == 0),vSNames]
    XR <- model.matrix(formula(paste("~", fixed.part)),
                       regR)
    R <- diag(sigma2R, nrow = nrow(regS), ncol = nrow(regS))
    G <- sigma(mEst) ^ 2 * crossprod(getME(mEst, "Lambdat"))
    div = eval(parse(text = paste(paste(
      "reg$", paste(division),
      sep = ""
    ))))
    Ysim = matrix(NA, nrow = N, ncol = L)
    Ysim[(con == 1),] = YS
    divS = div[con == 1]
    div = as.factor(div)
    divS = as.factor(divS)
    con0 <- (con == 0)
    D = nlevels(div)
    i = 1
    while (i <= D) {
      con_divS <- (divS == levels(div)[i])
      con_div <- (div == levels(div)[i])
      nd <- sum(con_divS)
      Nd <- sum(con_div)
      if (Nd > nd) {
        if (nd == 0) {
          Ysim[con_div,] <- Xbeta[con_div,] + replicate(L,
                                                        c(matrix(
                                                          rep(rnorm(1, 0, Sv), Nd), nrow = 1,
                                                          byrow = T
                                                        ))) + matrix(rnorm(L * Nd, 0, sigmaR),
                                                                     ncol = L)
        }
        else {
          Ysim[(con_div & con0),] <- matrix(rep((
            as.matrix(Xbeta[(con_div & con0),]) + matrix(S2v, Nd - nd, nd) %*% 
              solve(matrix(S2v, nd, nd) + sigma2R * diag(nd)) %*% 
              as.matrix(YS[con_divS] - Xbeta[(con_div &
                                                                                                                                           (con == 1)),])
          ), L), ncol = L) +
            replicate(L, c(matrix(
              rep(rnorm(1, 0, Sv),
                  (Nd - nd)), nrow = 1, byrow = T
            ))) + matrix(rnorm(L *
                                 (Nd - nd), 0, sigmaR), ncol = L)
        }
      }
      i = i + 1
    }
    
    if (is.null(backTrans)) {
      YbackTranssim <- Ysim
    }  else {
      YbackTranssim <- backTrans(Ysim)
    }
    
    thetaP <- sapply(1:L, function(i) {
      thetaFun(YbackTranssim[, i])
    })
    
    thetaP <- rowMeans(matrix(thetaP, ncol = L))
    
    outl <- list(
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
      weights = weights,
      Z = Z,
      ZBlockNames =  ZBlockNames,
      X = X,
      ZS = ZS,
      XR = XR,
      ZR = ZR,
      eS = eS,
      vS = vSDF,
      fixed.part = fixed.part,
      random.part = random.part,
      division = division,
      backTrans = backTrans,
      thetaFun = thetaFun,
      L = L
    )
    class(outl) = "ebpLMMne"
    return(outl)
  }
