lwzl <-
function(listRanef, reg) {
  llwzl <- lapply(1:length(listRanef), function(i){
    
    x <- listRanef[[i]]
    ranefname <- names(listRanef[i])
    k <- length(unique(reg[ ,ranefname])) #all v
    
    df <- data.frame(nrw = 1:nrow(x), rnames = rownames(x), x)
    nrwsamp <- as.vector(sample(df$nrw, k, replace = TRUE))
    
    dfsamp <- NULL
    for (i in 1:k){
      dfsamp <- rbind(dfsamp, df[nrwsamp[i], ])
    }
    
    raneftotal <- data.frame(
      ranef = as.vector(t(as.matrix(dfsamp[ ,-c(1,2)])))
    )
    return(list(raneftotal = raneftotal, ranefname = ranefname, k = k, df = df, 
                dfsamp = dfsamp))
  })
  names(llwzl) <- names(listRanef)
  
  tablwzl <-
    dplyr::bind_rows(lapply(1:length(listRanef), function(j) {
      data.frame(llwzl[[j]]$raneftotal)
    }))
  
  return(list(tablwzl = tablwzl, llwzl = llwzl))
}
