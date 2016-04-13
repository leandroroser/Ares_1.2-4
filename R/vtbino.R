
vtbino <- function(dens, Q=NULL, control=list()) {
  
  ctrl <- list(geps=1e-3,gmax=100,length=1000,span=0.05,eeps=1e-06,emax=5000)
  ctrl[names(control)]=control
  size <- length(dens)
  if(is.null(Q)) {
    Q <- emtbino(dens=dens,G=floor((size-1)/2),control=ctrl)
    Q <- mtbino(Q,span=ctrl$span)
    Q <- emtbino(dens=dens,Q=Q,control=ctrl)
  }
  
  step <- 1
  grid <- range(which(dens>0))
  grid <- seq(grid[1],grid[2],length=ctrl$length)
  grid <- ntbino(mu=grid,size=size)
  grad <- gtbino(dens=dens,Q=Q,prob=grid)
  i <- which.max(grad)
  while((grad[i]>ctrl$geps)&&(step<=ctrl$gmax)) {
    Q <- sltbino(dens=dens,Q=Q,prob=grid[i])
    Q <- emtbino(dens=dens,Q=Q,control=ctrl)
    grad <- gtbino(dens=dens,Q=Q,prob=grid)
    i <- which.max(grad)
    step <- step+1
  }
  
  Q <- mtbino(Q,span=0)
  attr(Q,"convergence")=(grad[i]<ctrl$geps)
  Q
}