EmpCM <-
function(model) {
  x <- ranef(model)
  L = length(x)
  centered <- lapply(1:L, function(i) {
    x[[i]] - colMeans(x[[i]])
  })
  empirical_covariance_matrix = lapply(1:L, function(i) {
    crossprod(as.matrix(centered[[i]])) / nrow(centered[[i]])
  })
  empirical_covariance_matrix
}
