
batbino <- function(dens,Q,total=NULL,control=list()) {
  
  ctrl <- list(eeps=1e-06,emax=5000)
  ctrl[names(control)]=control
  k <- nrow(Q)
  Q <- list(Q)
  if(k>1) {
    for (i in 2:k) {
      Q[[i]] <- mpair(Q=Q[[(i-1)]])
      Q[[i]] <- emtbino(dens=dens,Q=Q[[i]],control=ctrl)
    }
    Q <- rev(Q)
  }
  
  if(!is.null(total)) {
    size <- length(dens)
    attr(Q,"loglik") <- total*sapply(Q,function(D) sum(dens*log(dmtbino(x=1:size,size=size,Q=D))))
  }
  
  Q
}