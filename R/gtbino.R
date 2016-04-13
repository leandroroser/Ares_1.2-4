
gtbino <- function(dens,Q,prob) {
  size <- length(dens)
  i <- (1:size)[dens>0]
  prob <- sapply(i,dtbino,size=size,prob=prob)
  prob <- prob%*%(dens[i]/dmtbino(x=i,size=size,Q=Q))
  drop(prob)-1
}