summary.EBLUP <- 
  function(object, ...){
  cat(paste('Value of the predictor of the defined linear combination of the dependent variable =', round(object$thetaP, 4), 
             '\nto see the details, please use str()', '\n', '\n'))
    cat(paste('Sample size = ', length(object$YS), '\n'))
    cat(paste('Dataset size = ', nrow(object$reg), '\n', '\n'))
    print(summary(object$mEst))
  }
