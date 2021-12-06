###Running Means####
runningmeans <- function(x, start){
    warning("This function will be removed in version 0.9.1")
  n.iter=length(x)
  val=sapply(start:n.iter, function(ii)  mean(x[start:ii]))
  return(val)
}
