doubleBootFuture <-
function (predictor, B1, B2, p, q) 
{
  if (inherits(predictor, 'EBLUP') == F & inherits(predictor, 'plugInLMM') == F & 
    inherits(predictor, 'ebpLMMne') == F) {
    stop("wrong predictor")
  }
  if (B1 < 1) {
    stop("B1 must be > 0")
  }
  if (B2 < 1) {
    stop("B2 must be > 0")
  }
  B1 <- B1 + 1
  bootPar1 <- bootParFuture(predictor, B1, p)
  
  positiveDefiniteEstGlev1 <- bootPar1$positiveDefiniteEstG 
  
  if (positiveDefiniteEstGlev1 == FALSE) {
    cat(paste("non-positive definite estimated covariance matrix of random effects used at the first level - y is generated at the first level based on a model without random effects", "\n"))
  }
  
   
  number_of_predictors <- nrow(bootPar1$error)
  Ysim1 <- bootPar1$Ysim
  Ysim1S <- Ysim1[predictor$con == T, ]
  B1 <- B1 - 1
  error1 <- matrix(bootPar1$error[, 1:B1], nrow = number_of_predictors)
  estQAPE <- sapply(1:number_of_predictors, function(i) quantileNaN(abs(error1[i, 
  ]), probs = p))
  class <- class(predictor)
  
  if (class == "plugInLMM") {
    bootPar2 <- lapply(1:B1, function(i) {
      YS1 <- Ysim1S[, i]
      pred2 <- plugInLMM(YS1, predictor$fixed.part, predictor$random.part, 
                         predictor$reg, predictor$con, predictor$weights, 
                         predictor$backTrans, predictor$thetaFun)
      bootPar2 <- bootParFuture(pred2, B2, p)
      positiveDefiniteEstGlev2 <- bootPar2$positiveDefiniteEstG
      error2 <- bootPar2$error
      return(list(error2 = error2, positiveDefiniteEstGlev2 = positiveDefiniteEstGlev2))
    })
  }
  if (class == "ebpLMMne") {
    bootPar2 <- lapply(1:B1, function(i) {
      YS1 <- Ysim1S[, i]
      pred2 <- ebpLMMne(YS1, predictor$fixed.part, predictor$division, 
                        predictor$reg, predictor$con, predictor$backTrans, 
                        predictor$thetaFun, predictor$L)
      bootPar2 <- bootParFuture(pred2, B2, p)
      positiveDefiniteEstGlev2 <- bootPar2$positiveDefiniteEstG
      error2 <- bootPar2$error
      return(list(error2 = error2, positiveDefiniteEstGlev2 = positiveDefiniteEstGlev2))
    })
  }
  if (class == "EBLUP") {
    bootPar2 <- lapply(1:B1, function(i) {
      YS1 <- Ysim1[, i][predictor$con == T]
      pred2 <- EBLUP(YS1, predictor$fixed.part, predictor$random.part, 
                     predictor$reg, predictor$con, predictor$gamma, 
                     predictor$weights, estMSE = FALSE)
      bootPar2 <- bootParFuture(pred2, B2, p)
      positiveDefiniteEstGlev2 <- bootPar2$positiveDefiniteEstG
      error2 <- bootPar2$error
      return(list(error2 = error2, positiveDefiniteEstGlev2 = positiveDefiniteEstGlev2))
    })
  }
  
  
  positiveDefiniteEstGlev2 <- sum(as.numeric(sapply(1:B1, function(i) (bootPar2[[i]]$positiveDefiniteEstGlev2))))
  if (positiveDefiniteEstGlev2 < B1) {
    cat(paste("number of cases with non-positive definite estimated covariance matrix of random effects used at the second level:", (B1-positiveDefiniteEstGlev2), "out of", B1, "\n"))
  }
  
  error2 <- lapply(1:number_of_predictors, function(j) { 
    t(sapply(bootPar2, function(x) x$error2[j, ]))
  })
  
  mse4 <- sapply(1:number_of_predictors, function(i) sum(error1[i, 
  ]^2)/B1)
  mse7 <- sapply(1:length(error2), function(i) sum(error2[[i]]^2)/(B1 * 
                                                                     B2))
  mse8 <- 2 * mse4 - mse7
  u1_2 <- sapply(1:number_of_predictors, function(j) {
    sapply(1:B1, function(i) {
      u <- 2 * (error1[j, i]^2) - mean(error2[[j]][i, ]^2)
    })
  })
  u2_2 <- sapply(1:number_of_predictors, function(j) {
    sapply(1:B1, function(i) {
      u <- 2 * (error1[j, i]^2) - error2[[j]][i, 1]^2
    })
  })
  mse10 <- colMeans(matrix(u2_2, ncol = number_of_predictors))
  u3_2 <- sapply(1:number_of_predictors, function(j) {
    sapply(1:B1, function(i) {
      u <- error1[j, i]^2 + bootPar1$error[j, (i + 1)]^2 - 
        error2[[j]][i, 1]^2
    })
  })
  mse12 <- colMeans(matrix(u3_2, ncol = number_of_predictors))
  u1_2mod <- sapply(1:number_of_predictors, function(j) {
    sapply(1:B1, function(i) {
      u <- 2 * (error1[j, i]^2) - mean(error2[[j]][i, ]^2)
      ifelse(u < 0, error1[j, i]^2, u)
    })
  })
  u1_2mod <- matrix(u1_2mod, ncol = number_of_predictors)
  mse8mod <- colMeans(u1_2mod)
  u2_2mod <- sapply(1:number_of_predictors, function(j) {
    sapply(1:B1, function(i) {
      u <- 2 * (error1[j, i]^2) - error2[[j]][i, 1]^2
      ifelse(u < 0, error1[j, i]^2, u)
    })
  })
  u2_2mod <- matrix(u2_2mod, ncol = number_of_predictors)
  mse10mod <- colMeans(u2_2mod)
  u3_2mod <- sapply(1:number_of_predictors, function(j) {
    sapply(1:B1, function(i) {
      u <- error1[j, i]^2 + bootPar1$error[j, (i + 1)]^2 - 
        error2[[j]][i, 1]^2
      ifelse(u < 0, error1[j, i]^2, u)
    })
  })
  u3_2mod <- matrix(u3_2mod, ncol = number_of_predictors)
  mse12mod <- colMeans(u3_2mod)
  qape17 <- sapply(1:number_of_predictors, function(i) quantileNaN(sqrt(u1_2mod)[, 
                                                                              i], probs = p))
  qape18 <- sapply(1:number_of_predictors, function(i) quantileNaN(sqrt(u2_2mod)[, 
                                                                              i], probs = p))
  qape19 <- sapply(1:number_of_predictors, function(i) quantileNaN(sqrt(u3_2mod)[, 
                                                                              i], probs = p))
  mse14 <- ifelse(mse4 >= mse7, mse8, mse4 * exp((mse4 - mse7)/mse7))
  con1516 <- sapply(1:number_of_predictors, function(j) {
    (1/mse4[j]) * mean(error2[[j]][, 1]^2)
  })
  mse15 <- ifelse(con1516 < q, q * mse4, mse10)
  mse16 <- ifelse(con1516 < q, q * mse4, mse12)
  return(list(estMSE_param = mse4, estMSE_db_B2 = mse8, estMSE_db_B2_WDZ = mse8mod, 
              estMSE_db_B2_HM = mse14, estMSE_db_1 = mse10, estMSE_db_1_WDZ = mse10mod, 
              estMSE_db_1_EF = mse15, estMSE_db_telesc = mse12, estMSE_db_telesc_WDZ = mse12mod, 
              estMSE_db_telesc_EF = mse16, estQAPE_param = estQAPE, 
              estQAPE_db_B2 = qape17, estQAPE_db_1 = qape18, estQAPE_db_telesc = qape19, 
              error1 = error1, error2 = error2, corSquaredError1_db_B2 = u1_2, 
              corSquaredError1_db_1 = u2_2, corSquaredError1_db_telesc = u3_2, 
              corSquaredError1_db_B2_WDZ = u1_2mod, corSquaredError1_db_1_WDZ = u2_2mod, 
              corSquaredError1_db_telesc_WDZ = u3_2, positiveDefiniteEstGlev1 = positiveDefiniteEstGlev1, positiveDefiniteEstGlev2 = positiveDefiniteEstGlev2))
}
