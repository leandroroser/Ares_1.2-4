
mtbino <- function(Q,span) {
  Q <- Q[order(Q[,1]),,drop=FALSE]
  k <- nrow(Q)
  
  if(k>1) {
    l <- Q[-1,1]-Q[-k, 1]
    r <- (1:k)[c(l,2)>span]
    l <- (1:k)[c(2,l)>span]
    nQ <- matrix(0,length(l),2)
    for(i in seq(along=l)) {
      k <- l[i]:r[i]
      nQ[i,2] <- sum(Q[k,2])
      nQ[i,1] <- sum(Q[k,1]*Q[k,2])/nQ[i,2]
    }
    Q <- nQ
  }
  Q
}