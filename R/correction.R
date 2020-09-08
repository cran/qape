correction <-
function(model) { 
  R = EstCM(model)
  S = EmpCM(model)
  A = lapply(1:length(R), function(i) {
    LS = t(chol(S[[i]]))
    LR = t(chol(R[[i]]))
    t(LR %*% solve(LS))
  })
  A
}
