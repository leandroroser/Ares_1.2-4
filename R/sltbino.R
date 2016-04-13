
sltbino <- function(dens,Q,prob) {
  
  size <- length(dens)
  i <- which(dens>0)
  dens <- dens[i]
  r <- which.min(gtbino(dens=dens,Q=Q,prob=Q[,1]))
  f <- dmtbino(x=i,size=size,Q=Q)
  fadd <- dtbino(x=i,size=size,prob=prob)
  frem <- dtbino(x=i,size=size,prob=Q[r,1])
  a <- (fadd-frem)*Q[r,2]/f
  d1 <- c(sum(dens*a),sum(dens*a/(a+1)))
  d2 <- -c(sum(dens*a^2), sum(dens*(a/(a+1))^2))
  
  if(all(is.finite(c(d1,d2)))) {
    if((d2[1]<=d2[2])||(d2[1]<=(d1[2]-d1[1]))) {
      a=-d1[1]/d2[1]
    } else {
      a=-d1[1]/(d1[2]-d1[1])
    }	
    a=ifelse(a>1,NA,ifelse(a<0,NA,a))
  
  } else {
    a=NA
  }
  
  if(is.na(a)) {
    z=Q[r,2]*(fadd-frem)
    lik=function(p){ return(-sum(dens*log(f+p*z))) }
    a=suppressWarnings(optimize(lik,0:1)$minimum)
  }
  
  if(a>=1) {
    Q[r,1]=prob
  
  } else {
    Q=cbind(c(prob,Q[r,1],Q[-r,1]),c(a*Q[r,2],(1-a)*Q[r,2],Q[-r,2]))
  }
  
  Q
}