
dtbino <- function(x,size,prob) {
  z <- suppressWarnings(dbinom(x,size,prob)/(1-dbinom(0,size, prob)))
  ifelse(is.finite(z),z,ifelse(x==0,Inf,ifelse(x==1,1,0)))
}