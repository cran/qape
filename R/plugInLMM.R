plugInLMM <-
function(YS, fixed.part, random.part, reg, con, weights, backTrans, thetaFun){
  
  model <- formula(paste('YS', '~', fixed.part, '+', random.part))
  regS <- subset(reg, con == 1)
  regR <- subset(reg, con == 0)
  
  if (is.null(weights)) {
    weights <- rep(1, nrow(reg))
    weightS <- rep(1, nrow(regS))} else {
    weightS <- subset(weights, con == 1)}
  mEst <- lmer(model, weights = weightS, data.frame(YS, regS, weightS))
  tst <- unlist(apply((as.matrix(reg[ ,names(ranef(mEst))])),2,unique))
  if ((length(tst[duplicated(tst)])) > 0) {stop(paste("There are at least two 'random.part' variables with at least one the same value. 
                                                 Rename the values of 'random.part' variables - you can use qape::modifyDataset"))}
  Zobj <- Zfun(model, reg)
  Z <- Zobj$Z
  ZBlockNames <- Zobj$ZBlockNames
  vNames <- make.unique(Zfun(model, reg)$vNames,sep = ".") 
  colnames(Z) <- vNames
  X <- model.matrix(formula(paste('~', fixed.part)), reg)
  ZS <- getME(mEst, name = "Z")
  vSNames <- make.unique(colnames(ZS),sep = ".") 
  colnames(ZS) <- vSNames 
  ZR <- Z[(con == 0),vSNames]
  beta <- mEst@beta
  Xbeta <- X %*% beta
  XS <- getME(mEst, name = 'X')
  vS <- as.vector(getME(mEst, name = "b"))
  vSDF <- data.frame(vSNames, vS)
  eS <- residuals(mEst)
  XR <- model.matrix(formula(paste('~', fixed.part)), regR)
  sigma2R <- sigma(mEst)^2
  R <- diag(sigma2R/weightS, nrow = nrow(regS), ncol = nrow(regS))
  G <- sigma(mEst)^2 * crossprod(getME(mEst,"Lambdat")) 
  
  Y <- rep(NA,nrow(reg))
  Y[con == 1] <- YS
  YP <- Y[con == 0] <- XR %*% beta + ZR %*% vS
 
  if (is.null(backTrans))
    YbackTrans <- Y else
  YbackTrans <- backTrans(Y)

  if (is.null(backTrans))
    YPbackTrans <- YP else
      YPbackTrans <- backTrans(YP)

  thetaP <- thetaFun(YbackTrans)

  outl <- list(
    fixed.part = fixed.part,
    random.part = random.part, 
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
    weightS = weightS,
    Z = Z,
    ZBlockNames = ZBlockNames,
    X = X,
    ZS = ZS,
    XR = XR,
    ZR = ZR,
    eS = eS,
    vS = vSDF,
    thetaFun = thetaFun, 
    backTrans = backTrans)
    class(outl) = "plugInLMM"
  return(outl)
}
