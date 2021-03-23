Zfun <- function(model, data) 
{
  data <- data.frame(YY = rep(-1000, nrow(data)), data)
  names(data)[1] <- substring(model, 1)[2]
  obj <- mkReTrms(findbars(model), model.frame(subbars(model), data))
  Zt <- obj$Zt
  ZBlockNames <- make.unique(names(obj$cnms),sep = ".")
  vNames <- Zt@Dimnames[[1]]
  Zm <- as.matrix(Zt)
  Z <- t(Zm)
  vNames <- Zt@Dimnames[[1]]
  return(list(Z = as.matrix(Z), vNames = vNames, ZBlockNames = ZBlockNames))
}