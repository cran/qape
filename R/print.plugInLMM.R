print.plugInLMM <-
function (x, ...) {
  cat(paste(c('Value/s of the predictor of the defined function/s of the dependent variable =', 
              round(x$thetaP, 4)), collapse=" "))
  cat('\nto see the details, please use str()', '\n', '\n')
  cat(paste('Sample size = ', length(x$YS), '\n'))
  cat(paste('Dataset size = ', nrow(x$reg), '\n', '\n'))
}
