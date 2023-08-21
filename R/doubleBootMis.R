doubleBootMis <-
function(predictorLMM, predictorLMMmis, B1, B2, p, q) 
{
  
  if (inherits(predictorLMM, "plugInLMM") == F | inherits(predictorLMMmis, "plugInLMM") == F) {
    stop("wrong predictor/predictors")
  }
  
  if (!all.equal(predictorLMM$YS, predictorLMMmis$YS)){ 
    stop("'YS' object must be the same in both predictors")
  }
  
  if (!all.equal(predictorLMM$reg, predictorLMMmis$reg)){ 
    stop("'reg' object must be the same in both predictors")
  }
  
  if (!all.equal(predictorLMM$con, predictorLMMmis$con)){ 
    stop("'con' object must be the same in both predictors")
  }
  
  if (!all.equal(predictorLMM$thetaFun, predictorLMMmis$thetaFun)){ 
    stop("'thetaFun' object must be the same in both predictors")
  }
  
  if (B1 < 1) {
    stop("B1 must be > 0")
  }
  
  if (B2 < 1) {
    stop("B2 must be > 0")
  }
  
  B1 <- B1 + 1
  bootPar1 <- bootParMis(predictorLMM, predictorLMMmis, B1, p)

  positiveDefiniteEstGlev1 <- bootPar1$positiveDefiniteEstG 
  
  
  if (positiveDefiniteEstGlev1 == FALSE) {
        cat(paste("non-positive definite estimated covariance matrix of random effects used at the first level - y is generated at the first level based on a model without random effects", "\n"))
  }
  
   
  number_of_predictors <- nrow(bootPar1$errorLMM)
  Ysim1 <- bootPar1$Ysim
  Ysim1S <- Ysim1[predictorLMM$con == T, ]
  B1 <- B1 - 1
  error1LMM <- matrix(bootPar1$errorLMM[, 1:B1], nrow = number_of_predictors)
  error1LMMmis <- matrix(bootPar1$errorLMMmis[, 1:B1], nrow = number_of_predictors)
  estQAPElmm <- sapply(1:number_of_predictors, function(i) quantileNaN(abs(error1LMM[i, 
  ]), probs = p))
  estQAPElmmMis <- sapply(1:number_of_predictors, function(i) quantileNaN(abs(error1LMMmis[i, 
  ]), probs = p))
  
  bootPar2 <- lapply(1:B1, function(i) {
    YS1 <- Ysim1S[, i]
    pred2LMM <- plugInLMM(YS1, predictorLMM$fixed.part, predictorLMM$random.part, 
                          predictorLMM$reg, predictorLMM$con, predictorLMM$weights, 
                          predictorLMM$backTrans, predictorLMM$thetaFun)
    
    pred2LMMmis <- plugInLMM(YS1, predictorLMMmis$fixed.part, predictorLMMmis$random.part, 
                             predictorLMMmis$reg, predictorLMMmis$con, predictorLMMmis$weights, 
                             predictorLMMmis$backTrans, predictorLMMmis$thetaFun)
    
    bootPar2 <- bootParMis(pred2LMM, pred2LMMmis, B2, p)
    
    positiveDefiniteEstGlev2 <- bootPar2$positiveDefiniteEstG
    error2LMM <- bootPar2$errorLMM
    error2LMMmis <- bootPar2$errorLMMmis
    return(list(error2LMM = error2LMM, error2LMMmis = error2LMMmis, positiveDefiniteEstGlev2 = positiveDefiniteEstGlev2))
  })
  
 
  positiveDefiniteEstGlev2 <- sum(as.numeric(sapply(1:B1, function(i) (bootPar2[[i]]$positiveDefiniteEstGlev2))))
  
  if (positiveDefiniteEstGlev2 < B1) {
    cat(paste("number of cases with non-positive definite estimated covariance matrix of random effects used at the second level:", (B1-positiveDefiniteEstGlev2), "out of", B1, "\n"))
  }
  
  error2LMM <- lapply(1:number_of_predictors, function(j) { 
    t(sapply(bootPar2, function(x) x$error2LMM[j, ]))
  })
  
  error2LMMmis <- lapply(1:number_of_predictors, function(j) {
    t(sapply(bootPar2, function(x) x$error2LMMmis[j, ]))
  })
  
  mse4LMM <- sapply(1:number_of_predictors, function(i) sum(error1LMM[i, 
  ]^2)/B1)
  
  mse4LMMmis <- sapply(1:number_of_predictors, function(i) sum(error1LMMmis[i, 
  ]^2)/B1)
  
  mse7LMM <- sapply(1:length(error2LMM), function(i) sum(error2LMM[[i]]^2)/(B1 * 
                                                                              B2))
  
  mse7LMMmis <- sapply(1:length(error2LMMmis), function(i) sum(error2LMMmis[[i]]^2)/(B1 * 
                                                                           B2))
  
  mse8LMM <- 2 * mse4LMM - mse7LMM
  mse8LMMmis <- 2 * mse4LMMmis - mse7LMMmis
  
  u1_2LMM <- sapply(1:number_of_predictors, function(j) {
    sapply(1:B1, function(i) {
      u <- 2 * (error1LMM[j, i]^2) - mean(matrix(error2LMM[[j]],ncol=B2)[i, ]^2)
    })
  })
  
  u1_2LMMmis <- sapply(1:number_of_predictors, function(j) {
    sapply(1:B1, function(i) {
      u <- 2 * (error1LMMmis[j, i]^2) - mean(matrix(error2LMMmis[[j]], ncol=B2)[i, ]^2)  #############################
    })
  })
  
  u2_2LMM <- sapply(1:number_of_predictors, function(j) {
    sapply(1:B1, function(i) {
      u <- 2 * (error1LMM[j, i]^2) - matrix(error2LMM[[j]],ncol=B2)[i, 1]^2
    })
  })
  
  u2_2LMMmis <- sapply(1:number_of_predictors, function(j) {
    sapply(1:B1, function(i) {
      u <- 2 * (error1LMMmis[j, i]^2) - matrix(error2LMMmis[[j]],ncol=B2)[i, 1]^2
    })
  })
  
  mse10LMM <- colMeans(matrix(u2_2LMM, ncol = number_of_predictors))
  mse10LMMmis <- colMeans(matrix(u2_2LMMmis, ncol = number_of_predictors))
  
  u3_2LMM <- sapply(1:number_of_predictors, function(j) {
    sapply(1:B1, function(i) {
      u <- error1LMM[j, i]^2 + bootPar1$errorLMM[j, (i + 1)]^2 - 
        matrix(error2LMM[[j]],ncol=B2)[i, 1]^2
    })
  })
  
  u3_2LMMmis <- sapply(1:number_of_predictors, function(j) {
    sapply(1:B1, function(i) {
      u <- error1LMMmis[j, i]^2 + bootPar1$errorLMMmis[j, (i + 1)]^2 - 
        matrix(error2LMMmis[[j]],ncol=B2)[i, 1]^2
    })
  })
  
  mse12LMM <- colMeans(matrix(u3_2LMM, ncol = number_of_predictors))
  mse12LMMmis <- colMeans(matrix(u3_2LMMmis, ncol = number_of_predictors))
  
  u1_2modLMM <- sapply(1:number_of_predictors, function(j) {
    sapply(1:B1, function(i) {
      u <- 2 * (error1LMM[j, i]^2) - mean(matrix(error2LMM[[j]],ncol=B2)[i, ]^2)
      ifelse(u < 0, error1LMM[j, i]^2, u)
    })
  })
  u1_2modLMM <- matrix(u1_2modLMM, ncol = number_of_predictors)
  
  
  
  u1_2modLMMmis <- sapply(1:number_of_predictors, function(j) {
    sapply(1:B1, function(i) {
      u <- 2 * (error1LMMmis[j, i]^2) - mean(matrix(error2LMMmis[[j]],ncol=B2)[i, ]^2)
      ifelse(u < 0, error1LMMmis[j, i]^2, u)
    })
  })
  u1_2modLMMmis <- matrix(u1_2modLMMmis, ncol = number_of_predictors)
  
  mse8modLMM <- colMeans(u1_2modLMM)
  mse8modLMMmis <- colMeans(u1_2modLMMmis)
  
  u2_2modLMM <- sapply(1:number_of_predictors, function(j) {
    sapply(1:B1, function(i) {
      u <- 2 * (error1LMM[j, i]^2) - matrix(error2LMM[[j]],ncol=B2)[i, 1]^2
      ifelse(u < 0, error1LMM[j, i]^2, u)
    })
  })
  u2_2modLMM <- matrix(u2_2modLMM, ncol = number_of_predictors)
  
  
  u2_2modLMMmis <- sapply(1:number_of_predictors, function(j) {
    sapply(1:B1, function(i) {
      u <- 2 * (error1LMMmis[j, i]^2) - matrix(error2LMMmis[[j]],ncol=B2)[i, 1]^2
      ifelse(u < 0, error1LMMmis[j, i]^2, u)
    })
  })
  u2_2modLMMmis <- matrix(u2_2modLMMmis, ncol = number_of_predictors)  
  
  mse10modLMM <- colMeans(u2_2modLMM)
  mse10modLMMmis <- colMeans(u2_2modLMMmis)
  
  u3_2modLMM <- sapply(1:number_of_predictors, function(j) {
    sapply(1:B1, function(i) {
      u <- error1LMM[j, i]^2 + bootPar1$errorLMM[j, (i + 1)]^2 - 
        matrix(error2LMM[[j]],ncol=B2)[i, 1]^2
      ifelse(u < 0, error1LMM[j, i]^2, u)
    })
  })
  u3_2modLMM <- matrix(u3_2modLMM, ncol = number_of_predictors)
  
  u3_2modLMMmis <- sapply(1:number_of_predictors, function(j) {
    sapply(1:B1, function(i) {
      u <- error1LMMmis[j, i]^2 + bootPar1$errorLMMmis[j, (i + 1)]^2 - 
        matrix(error2LMMmis[[j]],ncol=B2)[i, 1]^2
      ifelse(u < 0, error1LMMmis[j, i]^2, u)
    })
  })
  u3_2modLMMmis <- matrix(u3_2modLMMmis, ncol = number_of_predictors)
  
  mse12modLMM <- colMeans(u3_2modLMM)
  mse12modLMMmis <- colMeans(u3_2modLMMmis)
  
  
  qape17LMM <- sapply(1:number_of_predictors, function(i) quantileNaN(sqrt(u1_2modLMM)[, 
                                                                                    i], probs = p))
  qape17LMMmis <- sapply(1:number_of_predictors, function(i) quantileNaN(sqrt(u1_2modLMMmis)[, 
                                                                                  i], probs = p))
  
  qape18LMM <- sapply(1:number_of_predictors, function(i) quantileNaN(sqrt(u2_2modLMM)[, 
                                                                                    i], probs = p))
  
  qape18LMMmis <- sapply(1:number_of_predictors, function(i) quantileNaN(sqrt(u2_2modLMMmis)[, 
                                                                                  i], probs = p))
  
  qape19LMM <- sapply(1:number_of_predictors, function(i) quantileNaN(sqrt(u3_2modLMM)[, 
                                                                                    i], probs = p))
  
  qape19LMMmis <- sapply(1:number_of_predictors, function(i) quantileNaN(sqrt(u3_2modLMMmis)[, 
                                                                                  i], probs = p))
  
  mse14LMM <- ifelse(mse4LMM >= mse7LMM, mse8LMM, mse4LMM * exp((mse4LMM - mse7LMM)/mse7LMM))
  mse14LMMmis <- ifelse(mse4LMMmis >= mse7LMMmis, mse8LMMmis, mse4LMMmis * exp((mse4LMMmis - mse7LMMmis)/mse7LMMmis))
  
  con1516LMM <- sapply(1:number_of_predictors, function(j) {
    (1/mse4LMM[j]) * mean(error2LMM[[j]][, 1]^2)
  })
  
  con1516LMMmis <- sapply(1:number_of_predictors, function(j) {
    (1/mse4LMMmis[j]) * mean(error2LMMmis[[j]][, 1]^2)
  })
  
  mse15LMM <- ifelse(con1516LMM < q, q * mse4LMM, mse10LMM)
  mse15LMMmis <- ifelse(con1516LMMmis < q, q * mse4LMMmis, mse10LMMmis)
  
  mse16LMM <- ifelse(con1516LMM < q, q * mse4LMM, mse12LMM)
  mse16LMMmis <- ifelse(con1516LMMmis < q, q * mse4LMMmis, mse12LMMmis)
  
  
  return(list(estMSE_param_LMM = mse4LMM, estMSE_db_B2_LMM = mse8LMM, estMSE_db_B2_WDZ_LMM = mse8modLMM, 
              estMSE_db_B2_HM_LMM = mse14LMM, estMSE_db_1_LMM = mse10LMM, estMSE_db_1_WDZ_LMM = mse10modLMM, 
              estMSE_db_1_EF_LMM = mse15LMM, estMSE_db_telesc_LMM = mse12LMM, estMSE_db_telesc_WDZ_LMM = mse12modLMM, 
              estMSE_db_telesc_EF_LMM = mse16LMM, estQAPE_param_LMM = estQAPElmm, 
              estQAPE_db_B2_LMM = qape17LMM, estQAPE_db_1_LMM = qape18LMM, estQAPE_db_telesc_LMM = qape19LMM, 
              error1LMM = error1LMM, error2LMM = error2LMM, corSquaredError1_db_B2_LMM = u1_2LMM, 
              corSquaredError1_db_1_LMM = u2_2LMM, corSquaredError1_db_telesc_LMM = u3_2LMM, 
              corSquaredError1_db_B2_WDZ_LMM = u1_2modLMM, corSquaredError1_db_1_WDZ_LMM = u2_2modLMM, 
              corSquaredError1_db_telesc_WDZ_LMM = u3_2modLMM,
              estMSE_param_LMMmis = mse4LMMmis, estMSE_db_B2_LMMmis = mse8LMMmis, estMSE_db_B2_WDZ_LMMmis = mse8modLMMmis, 
              estMSE_db_B2_HM_LMMmis = mse14LMMmis, estMSE_db_1_LMMmis = mse10LMMmis, estMSE_db_1_WDZ_LMMmis = mse10modLMMmis, 
              estMSE_db_1_EF_LMMmis = mse15LMMmis, estMSE_db_telesc_LMMmis = mse12LMMmis, estMSE_db_telesc_WDZ_LMMmis = mse12modLMMmis, 
              estMSE_db_telesc_EF_LMMmis = mse16LMMmis, estQAPE_param_LMMmis = estQAPElmmMis, 
              estQAPE_db_B2_LMMmis = qape17LMMmis, estQAPE_db_1_LMMmis = qape18LMMmis, estQAPE_db_telesc_LMMmis = qape19LMMmis, 
              error1LMMmis = error1LMMmis, error2LMMmis = error2LMMmis, corSquaredError1_db_B2_LMMmis = u1_2LMMmis, 
              corSquaredError1_db_1_LMMmis = u2_2LMMmis, corSquaredError1_db_telesc_LMMmis = u3_2LMMmis, 
              corSquaredError1_db_B2_WDZ_LMMmis = u1_2modLMMmis, corSquaredError1_db_1_WDZ_LMMmis = u2_2modLMMmis, 
              corSquaredError1_db_telesc_WDZ_LMMmis = u3_2modLMMmis, positiveDefiniteEstGlev1 = positiveDefiniteEstGlev1, positiveDefiniteEstGlev2 = positiveDefiniteEstGlev2))
}
