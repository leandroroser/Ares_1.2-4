
omtbino <- function(size,Q) {
  sum(1/(1/(1-Q[,1])^size-1)*Q[,2])
}