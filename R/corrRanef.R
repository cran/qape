corrRanef <-
function(model) {
  obj <- correction(model)
    if (inherits(obj, "try-error") == F) {
  rnames <- names(ranef(model))
  L <- length(ranef(model))
  lranef <- lapply(1:L, function(i) {
    as.data.frame(as.matrix((ranef(model)[[i]])) %*%
                    as.matrix(correction(model)[[i]]))
  })
  names(lranef) <- rnames}
  else {
  cat(paste("the correction of random effects cannot be made, random effects without correction are computed", "\n"))
  lranef <-ranef(model)}
  return(lranef)
}
