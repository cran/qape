print.EBLUP <- 
  function(x, ...){
  cat(paste('Value/s of the predictor =', round(x$thetaP, 4), 
             '\nto see the details, please use str()', '\n', '\n'))
    cat(paste('Sample size = ', length(x$YS), '\n'))
    cat(paste('Dataset size = ', nrow(x$reg), '\n', '\n'))
     }
