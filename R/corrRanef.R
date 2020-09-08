corrRanef <-
function(model) {
  rnames <- names(ranef(model))
  
  L <- length(ranef(model))
  lranef <- lapply(1:L, function(i) {
    as.data.frame(as.matrix((ranef(model)[[i]])) %*%
                    as.matrix(correction(model)[[i]]))
  })
  names(lranef) <- rnames
  return(lranef)
}
