normCholTest <- function(model,normTest){
  normTest <- match.fun(normTest)
  y <- as.matrix(lme4::getME(model, name = "y"))
  Z <- as.matrix(lme4::getME(model, name = "Z"))
  X <- as.matrix(lme4::getME(model, name = "X"))
  eG <- crossprod(lme4::getME(model,"Lambdat")) * sigma(model)^2
  eR <- diag(sigma(model)^2, nrow = nrow(X), ncol = nrow(X))
  
  eV <- as.matrix(Z %*% eG %*% t(Z) + eR)
  inveV <- solve(eV) #inverse of estimated V
  
  eBeta <- as.matrix(lme4::getME(model, name = "beta"))
  res <- y - X %*% as.matrix(eBeta)
  transformed_residuals <- t(chol(inveV)) %*% res
  testResults <- normTest(transformed_residuals)
  return(list(test = testResults, p.value = testResults$p.value))
         }
