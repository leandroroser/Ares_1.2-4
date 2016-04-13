

rpldtbino <- function(Q,size,to.conditional=FALSE) {
  x <- 1-(1-Q[,1])^size
  
  if(to.conditional) {
    Q[,2] <- Q[,2]*x
  
  } else {
    Q[,2] <- Q[,2]/x
    Q[,2] <- Q[,2]/sum(Q[,2])
  }
  
  Q
}