EBLUP <- function(YS, fixed.part, random.part, reg, con, gamma, 
                   weights, estMSE) {
  model <- formula(paste("YS", "~", fixed.part, 
                         "+", random.part))
  regS <- subset(reg, con == 1)
  regR <- subset(reg, con == 0)
  gammaS <- subset(gamma, con == 1)
  gammaR <- subset(gamma, con == 0)
  if (is.null(weights)) {
    weights <- rep(1, nrow(reg))
    weightS <- rep(1, nrow(regS))}  else {
    weightS <- subset(weights, con == 1)}
  mEst <- lmer(model, weights = weightS, data.frame(YS, regS, weightS))
  tst <- unlist(apply(as.matrix(reg[,names(ranef(mEst))]),2,unique))
  if ((length(tst[duplicated(tst)])) > 0) {stop(paste("There are at least two 'random.part' variables with at least one the same value. 
                                                 Rename the values of 'random.part' variables - you can use qape::modifyDataset"))}
  Zobj <- Zfun(model, reg)
  Z <- Zobj$Z
  ZBlockNames <- Zobj$ZBlockNames
  vNames <- make.unique(Zfun(model, reg)$vNames,sep = ".") 
  colnames(Z) <- vNames
  X <- model.matrix(formula(paste("~", fixed.part)), reg)
 
  ZS <- getME(mEst, name = "Z")
  vSNames <- make.unique(colnames(ZS),sep = ".") 
  colnames(ZS) <- vSNames 
 
  ZR <- Z[(con == 0),vSNames]
  beta <- mEst@beta
  Xbeta <- X %*% beta
  XS <- getME(mEst, name = "X")
  vS <- as.vector(getME(mEst, name = "b"))
  vSDF <- data.frame(vSNames, vS)
  eS <- residuals(mEst)
  XR <- model.matrix(formula(paste("~", fixed.part)), 
                     regR)
  thetaP <- as.numeric(gammaS %*% YS + gammaR %*% XR %*% beta + 
                         gammaR %*% ZR %*% as.vector(vS))
  sigma2R <- sigma(mEst)^2
  R <- diag(sigma2R/weightS, nrow = nrow(regS), ncol = nrow(regS))

  G <- sigma(mEst)^2 * crossprod(as.matrix(getME(mEst, "Lambdat")))
  if (estMSE == TRUE) {
  weightR <- subset(weights, con == 0)
  Rrr <- diag(sigma2R/weightR, nrow = nrow(regR), ncol = nrow(regR))
  invVss <- solve(ZS %*% G %*% t(ZS) + R)
  Vrr <- ZR %*% G %*% t(ZR) + Rrr
  Vrs <- ZR %*% G %*% t(ZS)
  VrsinvVss <- Vrs %*%  invVss
  VrsinvVssXs <- VrsinvVss %*% XS 
  g1 <- as.numeric(t(gammaR) %*% (Vrr - VrsinvVss %*% t(Vrs)) %*% gammaR)
  g2a <- t(gammaR) %*% (XR - VrsinvVssXs)
  g2 <- as.numeric(g2a %*% solve(t(XS) %*% invVss %*% XS) %*% t(g2a))
  neMSE <- g1 + g2
  } else{g1 <- g2 <- neMSE <- "not computed"}
  outl <- list(fixed.part = fixed.part, random.part = random.part, 
               thetaP = thetaP, beta = beta, Xbeta = Xbeta, sigma2R = sigma2R, 
               R = R, G = G, model = model, mEst = mEst, YS = YS, reg = reg, 
               con = con, regS = regS, regR = regR, gamma = gamma, 
               gammaS = gammaS, gammaR = gammaR, weights = weights, Z = Z, 
               ZBlockNames = ZBlockNames, X = X, ZS = ZS, 
               XR = XR, ZR = ZR, eS = eS, vS = vSDF, g1 = g1, g2 = g2, 
               neMSE = neMSE)
  class(outl) = "EBLUP"
  return(outl)
}