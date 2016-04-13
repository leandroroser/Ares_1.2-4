
dmtbino <- function(x,size,Q) {
  drop(sapply(Q[,1],dtbino,x=x,size=size)%*%Q[,2])
}