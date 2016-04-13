

fmtbino <- function(count=NULL,x=NULL,size=NULL) {
  if(is.null(count)) {
    x <- as.data.frame(table(x[x>0]))
    x[,1] <- unclass(x[,1])
    count <- rep(0,size)
    count[x[,1]] <- x[,2]
  }
  
  total <- sum(count)
  attr(count,"total") <- total
  
  count/total
}