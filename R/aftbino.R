
aftbino <- function(size,total,Q,sizes=NULL) {
  
  Q <- mtbino(Q,span=0)
  
  if(is.null(sizes)) {
    sizes <- 1:(4*size)
  }
  
  if(Q[1,1]>0) {
    out <- Q[,2]/(1-(1-Q[,1])^size)
    out <- drop(out%*%(1-outer(1-Q[,1],sizes,"^")))
  
  } else {
    out <- Q[-1,2]/(1-(1-Q[-1,1])^size)
    out <- drop(out%*%(1-outer(1-Q[-1,1],sizes,"^")))
    out <- out+Q[1,2]*sizes/size
  }
  
  out <- cbind(size=sizes,esti=total*out)
  
  invisible(out)
}