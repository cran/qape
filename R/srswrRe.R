srswrRe <- function(listRanef, reg) 
{
  lsrswrRe <- lapply(1:length(listRanef), function(i) {
    x <- listRanef[[i]]
    ranefname <- names(listRanef[i])
    k <- length(unique(reg[, ranefname]))
    df <- data.frame(nrw = 1:nrow(x), rnames = rownames(x), 
                     x)
    nrwsamp <- as.vector(sample(df$nrw, k, replace = TRUE))
    dfsamp <- NULL
    for (i in 1:k) {
      dfsamp <- rbind(dfsamp, df[nrwsamp[i], ])
    }
    raneftotal <- data.frame(ranef = as.vector(t(as.matrix(dfsamp[, 
                                                                  -c(1, 2)]))))
    return(list(raneftotal = raneftotal, ranefname = ranefname, 
                k = k, df = df, dfsamp = dfsamp))
  })
  names(lsrswrRe) <- names(listRanef)
  tablsrswrRe <- dplyr::bind_rows(lapply(1:length(listRanef), 
          function(j) {
          reNames <- unique(reg[, names(listRanef[j])])
          numberOfLevels <- length(reNames)
          tim <- length(unlist(lsrswrRe[[j]]$raneftotal)) / numberOfLevels
          refNames <- make.unique(as.character(matrix(rep(reNames,tim),
                                          nrow = tim, byrow = TRUE)), sep = ".")
          data.frame(refNames,lsrswrRe[[j]]$raneftotal)
                                         }))
  return(list(tablsrswrRe = tablsrswrRe, lsrswrRe = lsrswrRe))
}