
removetsp <- function(Q,epsilon=1e-6) {
  Q <- Q[Q[,2]>epsilon,,drop=FALSE]
  Q[,2] <- Q[,2]/sum(Q[,2])
  Q
}