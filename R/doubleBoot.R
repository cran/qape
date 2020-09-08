doubleBoot <- function(predictor, B1, B2, p, q){ 
  
  if (class(predictor) != 'EBLUP' & class(predictor) != 'plugInLMM'& class(predictor) != 'ebpLMMne') { ### TZ zmieni?em
    stop("wrong predictor")
  }
  
  if (B1 < 1) { 
    stop("B1 must be > 0")
  }
  
  if (B2 < 1) { 
    stop("B2 must be > 0")
  }
  

##first level
B1 <- B1+1  
bootPar1 <- bootPar(predictor, B1, p)
number_of_predictors=nrow(bootPar1$error)

Ysim1 <- bootPar1$Ysim 
Ysim1S=Ysim1[predictor$con == T,] 

B1 <- B1-1
error1 <- matrix(bootPar1$error[,1:B1],nrow=number_of_predictors)

estQAPE = sapply(1:number_of_predictors, function(i) quantile(abs(error1[i,]), probs = p))

class <- class(predictor)


if (class == 'plugInLMM') {
  ##second level
  bootPar2 <- lapply(1:B1, function(i){
    YS1 <- Ysim1S[ ,i] 
    pred2 <- plugInLMM(YS1, predictor$fixed.part, predictor$random.part, predictor$reg, predictor$con, 
                       predictor$weights, predictor$backTrans, predictor$thetaFun)
    bootPar2 <- bootPar(pred2, B2, p)
    error2 <-  bootPar2$error 
    return(error2)
  })
}

if (class == 'ebpLMMne') {
  ##second level
  bootPar2 <- lapply(1:B1, function(i){
    YS1 <- Ysim1S[ ,i] 
    pred2 <- ebpLMMne(YS1, predictor$fixed.part, predictor$division, predictor$reg, predictor$con, 
                      predictor$backTrans, predictor$thetaFun, predictor$L)
    bootPar2 <- bootPar(pred2, B2, p)
    error2 <-  bootPar2$error 
    return(error2)
  })
}


if (class == 'EBLUP') {
##second level
bootPar2 <- lapply(1:B1, function(i){
  YS1 <- Ysim1[ ,i][predictor$con == T]
  pred2 <- EBLUP(YS1, predictor$fixed.part, predictor$random.part, predictor$reg, predictor$con, predictor$gamma, predictor$weights)
  
  bootPar2 <- bootPar(pred2, B2, p)
  error2 <-  bootPar2$error 
  return(error2)
})
}

error2=lapply(1:number_of_predictors,function(j){
  t(sapply(bootPar2,function(x)x[j,]))
})

mse4 <- sapply(1:number_of_predictors, function(i) sum(error1[i,]^2)/B1)
  
mse7 <- sapply(1:length(error2),function(i) sum(error2[[i]]^2)/(B1*B2))

mse8 <- 2*mse4 - mse7

u1_2 <- sapply(1:number_of_predictors, function(j){
  sapply(1:B1, function(i){
  u <- 2*(error1[j,i]^2) - mean(error2[[j]][i,]^2)    
})})

colMeans(u1_2) 

u2_2 <- sapply(1:number_of_predictors, function(j){
  sapply(1:B1, function(i){ 
  u <- 2*(error1[j,i]^2) - error2[[j]][i,1]^2
})}) 

mse10 <- colMeans(u2_2)

u3_2 <- sapply(1:number_of_predictors, function(j){
  sapply(1:B1, function(i){
  u <- error1[j,i]^2 + bootPar1$error[j,(i+1)]^2  - error2[[j]][i,1]^2
})}) 

mse12 <- colMeans(u3_2)

u1_2mod <- sapply(1:number_of_predictors, function(j){
  sapply(1:B1, function(i){  ### TZ tu i w analogicznych miejscach poni?ej by?o do (B-1), ale na g?rze zdefiniowa?a? we wcze?niejszej wersji B=B-1. Teraz jest do B1 i to jest OK.
    u <- 2*(error1[j,i]^2) - mean(error2[[j]][i,]^2)    
    ifelse(u<0,error1[j,i]^2, u)    
  })})

mse8mod <- colMeans(u1_2mod) 

u2_2mod <- sapply(1:number_of_predictors, function(j){
  sapply(1:B1, function(i){ 
    u <- 2*(error1[j,i]^2) - error2[[j]][i,1]^2
    ifelse(u<0,error1[j,i]^2, u)  
  })}) 

mse10mod <- colMeans(u2_2mod) 


u3_2mod <- sapply(1:number_of_predictors, function(j){
  sapply(1:B1, function(i){
    u <- error1[j,i]^2 + bootPar1$error[j,(i+1)]^2  - error2[[j]][i,1]^2
    ifelse(u<0,error1[j,i]^2, u) 
  })}) 


mse12mod <- colMeans(u3_2mod) 

qape17 <- sapply(1:number_of_predictors, function(i) quantile(sqrt(u1_2mod)[,i], probs = p))
qape18 <- sapply(1:number_of_predictors, function(i) quantile(sqrt(u2_2mod)[,i], probs = p))
qape19 <- sapply(1:number_of_predictors, function(i) quantile(sqrt(u3_2mod)[,i], probs = p))
mse14 <- ifelse(mse4 >= mse7, mse8, mse4*exp((mse4 - mse7)/mse7))

con1516 <- sapply(1:number_of_predictors, function(j) {
  (1/mse4[j])*mean(error2[[j]][,1]^2)
  })

mse15 <- ifelse(con1516 < q, q*mse4, mse10) 
mse16 <- ifelse(con1516 < q, q*mse4, mse12)

return(list(estMSE_param = mse4, estMSE_db_B2 = mse8, estMSE_db_B2_WDZ = mse8mod, estMSE_db_B2_HM = mse14,
            estMSE_db_1 = mse10, estMSE_db_1_WDZ = mse10mod, estMSE_db_1_EF = mse15,
            estMSE_db_telesc = mse12, estMSE_db_telesc_WDZ = mse12mod, estMSE_db_telesc_EF = mse16, 
            estQAPE_param=estQAPE, estQAPE_db_B2 = qape17, estQAPE_db_1  = qape18, estQAPE_db_telesc  = qape19
))


}

