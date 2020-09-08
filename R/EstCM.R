EstCM <-
function(model) {
    temp = VarCorr(model)
    if (length(temp) > 1){
      x = lapply(1:length(temp), function(i) {
        as.matrix(as.data.frame(VarCorr(model)[[i]]))
      }) } else {
        x = list(as.matrix(as.data.frame(VarCorr(model)[[1]])))
      }
    names(x) <- sub("\\..*", "", names(temp))
    
    lvar <- list()
    for (i in 1:length(ranef(model))){
      tmp <- list()
      
      cond <- names(ranef(model))[i]
      for (j in 1:length(x)){
        if (names(x)[j] == cond) tmp[j] <- x[j]
      }
      tmp <- plyr::compact(tmp)
      
      lvar[i] <- bdiag(tmp)
      names(lvar)[i] <- cond
      
    }
    
  return(estimated_covariance_matrix = lvar)
  }