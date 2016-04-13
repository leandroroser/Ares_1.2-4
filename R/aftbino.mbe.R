
aftbino.mbe <- function(count,estimated.richness=NULL,conf.level=0.95) {
  
  size <- length(count)
  if (is.null(estimated.richness)) {
    estimated.richness <- ifelse(all(count[1:2]>0), 
                                 count[1]^2/count[2]/2/size*(size - 1)+sum(count), 
                                 Inf)
  }
  posi <- (1:size)[count>0]
  count <- count[count>0]
  mat <- 1-outer(size-(1:size),posi, choose)/outer(rep(size,size),posi,choose)
  esti <- drop(mat%*%count)
  se <- drop(mat^2%*%count)-esti^2/estimated.richness
  se <- suppressWarnings(sqrt(se))*qnorm((1+conf.level)/2)
  out <- cbind(size=1:size,esti=esti,lwr=esti-se, upr=esti+se)
  
  invisible(out)
}
