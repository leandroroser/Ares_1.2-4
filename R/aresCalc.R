# Allelic richness estimation, with extrapolation beyond the sample size

aresCalc <- function(aresdata, pop = NULL, bootsize=NULL, maxsize=NULL) {
  
  if(class(aresdata) == "ARES") {
  if(is.null(pop)) {
      stop("please provide a pop argument")
  }
  rawdata <- attr(aresdata@OUT[[pop]], "output")
  
  } else if(class(aresdata) == "matrix" | class(aresdata) == "data.frame") {
    rawdata <- aresdata
  } else {
    stop("non valid class for aresdata")
  }
    
  alri.est <- 0
  
  # bootsize = bootstrap size
  if(is.null(bootsize))
  {
    bootsize <- 200
  }
  
  # how far to extrapolate
  if(is.null(maxsize))
  {
    maxsize <- 100
  }
  
  # convert the data to binary matrix
  rawdata  <-  rawdata > 0
  
  # convert the data to 'dens', which is used by the other functions
  dens <- fmtbino(x=rowSums(rawdata),size=ncol(rawdata))
  size <- length(dens)
  
  # Compute the 'allelic accumulation curve' from one 
  # to the nr of observed species, with moment-based estimates.
  total <- attr(dens, "total")
  est.mbe <- aftbino.mbe(count = dens * total)
  
  if(size >= maxsize) {
    alri.tot <- est.mbe[1:maxsize,2:4]
  
  } else {
    # The calculations to extrapolate & get confidence bounds
    Q <- vtbino(dens)
    
    # remove zeroprob from Q, to prevent an infinity returned by omtbino
    Qp <- Q
    Qp[,1] <- pmax(1e-8,Qp[,1])
    richness <- round((omtbino(size, Qp) + 1) * total)
    
    # extra safety, but probably a richness of one million does not make sense
    # anyway - to be chaned in a future version.
    if (is.infinite(richness)|richness>1e6) {
      richness <- 1e6
    }
    
    nvec <- rbinom(bootsize, richness, total/richness)
    prob <- dmtbino(1:size, size, Q)
    dmat <- sapply(nvec, function(x) rmultinom(1, x, prob))
    dmat <- dmat/matrix(nvec, size, bootsize, byrow = TRUE)
    Qlist <- apply(dmat, 2, vtbino)
    while (length(Qlist) < bootsize) {
      Qlist <- apply(dmat, 2, vtbino)
    }
    estimat <- matrix(0, maxsize, bootsize)
    
    for (i in 1:bootsize) {
      estimat[, i] <- aftbino(size, total = nvec[i], Q = Qlist[[i]], sizes = 1:maxsize)[,2]
    }

    ci <- t(apply(estimat, 1, quantile, c(0.025, 0.975)))
    esti <- aftbino(size, total = total, Q = Q, sizes = 1:maxsize)
    est.mle <- cbind(esti[, 2], ci)
    # glue the moment-based and maximum-likelyhood estimation together
    alri.tot <- rbind(est.mbe[,2:4], est.mle[(size+1):maxsize,1:3])
  }

  alri.tot
}