modifyDataset <- function(data, names) {
  tst <- unlist(apply(as.matrix(data[,names]),2,unique))
  if ((length(tst[duplicated(tst)])) == 0) {stop(paste("Values of the variables 
     listed in 'names' are unique - the dataset has not been modified."))}
  for (i in 1:length(names)) {
    prefix <- names[i]
    data[,names[i]] <- paste(prefix,data[,names[i]],sep = ".")
  }
  data  
}
