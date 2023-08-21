quantileNaN <- function (x, probs) {
  if (sum(is.nan(x)) > 0) rep(NaN,length(probs)) else {quantile(x, probs)}}
