Zfun <-
function(model, data){
        data <- data.frame(YY = rep(-1000, nrow(data)), data)
        names(data)[1] <- substring(model, 1)[2]
        Zt <- mkReTrms(findbars(model), model.frame(subbars(model), data))$Zt
        vNames <- Zt@Dimnames[[1]]
        Zm <- as.matrix(Zt)
        Z <- t(Zm)
        vNames <- Zt@Dimnames[[1]]
        return(list(Z = as.matrix(Z), vNames = vNames))
    }
