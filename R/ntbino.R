#' ntbino
#' @keywords internal

ntbino <- function(mu, size, neps=1e-8) {
  out <- rep(0,length(mu))
  out[mu>=size]=1
  i <- (mu>1)&(mu<size)
  
  if(any(i)) {
    mu <- mu[i]
    p <- mu/size
    d <- rep(2, length(p))
    
    while(max(abs(d))>neps) {
      d <- mu*(1-p)^(size-1)
      d <- (size*p-mu+(1-p)*d)/size/(1-d)
      p <- p-d
    }	
    
    out[i] <- p
  }
  
  out
}