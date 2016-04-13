

emtbino <- function(dens,G=1,Q=NULL,control=list()) {
  
  ctrl <- list(eeps=1e-06,emax=5000)
  ctrl[names(control)]=control
  size <- length(dens)
  i <- which(dens>0)
  dens <- dens[i]
  
  if(is.null(Q)) {
    if(G>1) {
      Q <- ntbino(mu=range(i),size=size)
      Q <- cbind(sort(runif(G,Q[1],Q[2])),(G:1)/sum(1:G))
    } else {
      return(cbind(ntbino(mu=sum(dens*i),size=size),1))
    }
  }
  
  like <- sum(dens*log(dmtbino(x=i,size=size,Q=Q)))
  
  for(j in 1:ctrl$emax) {
    mat <- rbind(sapply(i,dtbino,prob=Q[,1],size=size))
    mat <- t(mat*Q[,2])
    mat <- mat*(dens/rowSums(mat))
    Q[,2] <- colSums(mat)
    Q[,1] <- ntbino(drop(i%*%mat)/Q[,2],size)
    nlike <- sum(dens*log(dmtbino(x=i,size=size,Q=Q)))
    if(nlike-like<ctrl$eeps) {
      break
    } else {
      like = nlike
    }
  }
  
  Q
}