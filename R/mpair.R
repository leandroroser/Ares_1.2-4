
mpair <- function(Q) {
  Q <- Q[order(Q[,1]),,drop=FALSE]
  i <- nrow(Q)
  
  if(i>1) {
    i <- which.min(Q[-1,1]-Q[-i,1])
    Q[i,1] <- (Q[i,1]*Q[i,2]+Q[i+1,1]*Q[i+1,2])/(Q[i,2]+Q[i+1,2])
    Q[i,2] <- Q[i,2]+Q[i+1,2]
    Q <- Q[-(i + 1),,drop=FALSE]
  }
  
  Q
}