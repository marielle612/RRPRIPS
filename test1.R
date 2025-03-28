# S T A R T 

library(dr)
library(mbrdr)
library(ggplot2)
library(Renvlp)

##############################################
### Required Functions
##############################################
kmeans.min <- function(y, k=2, min=5){
  k.means <- kmeans(y,k) #k-means
  less.ct <- which(k.means$size<min) #less.ct: if cluster size < min(mininum size)
  if(length(less.ct)==0){
    new.cluster <- k.means$cluster
    new.size    <- table(k.means$cluster)
    n.slice     <- length(new.size)
  }else{
    repeat {
      s.cluster  <- k.means$cluster[!k.means$cluster %in% less.ct] # selection of cluster (size>=min)
      s.ncluster <- ifelse(length(table(s.cluster))<1,1, length(table(s.cluster)))# #of selectin(size>=min)
      new.kmeans <- kmeans(y, s.ncluster) #new k-means with s.ncluster
      
      new.cluster <- new.kmeans$cluster
      new.size    <- table(new.cluster)
      n.slice     <- length(new.size)
      
      less.rp <- which(new.kmeans$size<min) #less.rp : if new cluster size <min 
      if(length(less.rp)==0){ break }
      else{k.means <- kmeans(y, ifelse(n.slice-1<1, 1, n.slice-1))}
    }
  }
  myoutput        <- list(new.cluster, new.size, n.slice)
  names(myoutput) <- c("slice.indicator", "slice.sizes", "nslices")
  return(myoutput)
}

tr.vcc <- function (m1, m2) {
  m1 <- as.matrix(m1)
  m2 <- as.matrix(m2)
  k <- dim(m1)[2]
  P1 <-  m1%*% solve(t(m1) %*% m1)  %*% t(m1)
  P2 <-  m2%*% solve(t(m2) %*% m2)  %*% t(m2)
  
  P <- P1 %*% P2
  evals <- eigen(P)$values[1:k]
  r <- sqrt(sum(evals)/k)
  r <- 1-r
  q <- sqrt(prod(evals))
  q <- 1-q
  ans <- list(r=r, q=q)
  ## r and q: smaller better!!
  return(ans)
}
#######################################################################
#######################################################################

dr.wts <- function(object) {object$weights}

#######################################################################
#######################################################################
## Multiveraite version of Reduced-Rank Response Method
#######################################################################
#######################################################################

#######################################################################
## Multiveraite version of Reduced-Rank Response Combination (RRRComb)
#######################################################################
dr.M.rrrcomb<-function(object, nslices=2, dy=2, nclust=4, minc=5, slice.info=NULL,...) {
  z <- dr.z(object)
  x <- dr.x(object)
  y <- dr.y(object)
  wts <- dr.wts(object)
  
  r <- dim(y)[2]
  p <- dim(z)[2]
  
  ## Stating reduced-rank response mean method
  cand.mat <- NULL
  s0 <- mbrdr(y~x, method="prr")
  cand.mat <- cbind(cand.mat, s0$evectors[,1:dy])
  
  s1 <- update(s0, method="pfrr")
  cand.mat <- cbind(cand.mat, s1$evectors[,1:dy])
  
  s2 <- update(s0, method="upfrr")
  cand.mat <- cbind(cand.mat, s2$evectors[,1:dy])
  
  s3 <- update(s0, method="yc")
  cand.mat <- cbind(cand.mat, s3$evectors[,1:dy])
  
  pr.iter <- dim(cand.mat)[2]
  M<-0
  for (i in 1:pr.iter){ 
    t.i <- cand.mat[,i]
    t.i <- t.i /sqrt(sum(t.i^2))
    t.i.y <- c(y %*% t.i)
    ols1 <- dr.directions(dr( t.i.y~x, method="ols"))[,1]
    phd1 <- dr.directions(dr( t.i.y~x, method="phdres"))[,1]
    PI <- cbind(ols1, phd1)  
    
    if (nclust == 4) { kmeans.out <- dr.slices(PI, nslices=c(2,2))
    actual.k <- kmeans.out$nslices
    C.x <- kmeans.out$slice.indicator
    size.Cx <- kmeans.out$slice.sizes
    total.Cx <- sum(size.Cx) }  else { kmeans.out <- kmeans.min(PI, k=nclust, min=minc)
    actual.k <- kmeans.out$nslices
    C.x <- kmeans.out$slice.indicator
    size.Cx <- kmeans.out$slice.sizes
    total.Cx <- sum(size.Cx)   }
    
    M.j <- NULL
    for(k in 1:actual.k){
      sel.k <- C.x==k
      z.k <- z[sel.k, ] ;  y.k <- t.i.y[sel.k]
      wts.k <- wts[sel.k]
      
      slices <- dr.slices.arc(y.k, nslices)
      zmeans <- matrix(0,slices$nslices,NCOL(z.k))
      slice.weight <- rep(0,slices$nslices)  # NOT rep(0,NCOL(z.k))
      
      # compute weighted means within slice 
      wmean <- function (x, wts) { sum(x * wts) / sum (wts) }
      for (j in 1:slices$nslices){
        sel <- slices$slice.indicator==j
        zmeans[j,]<- apply(z.k[sel,],2, wmean, wts.k[sel])
        slice.weight[j]<-sum(wts.k[sel])}
      
      # get M matrix for sir
      ini.M <- t(apply(zmeans, 2, "*", sqrt(slice.weight))/ sqrt(sum(slice.weight)))
      M.j <- cbind(M.j, ( ini.M*(size.Cx[k]/total.Cx) ))
    } ## coordinate mean
    
    ##    if (dim(M.j)[2] != (nslices*nclust) ) {
    ##         diff.col <-  (nslices*nclust) - dim(M.j)[2]
    ##         M.j.zeros <- matrix( rep(0, (p*diff.col) ), c(p, diff.col) )
    ##         M.j <- cbind(M.j, M.j.zeros) }
    
    M.j <- M.j%*% t(M.j)
    M <- M+M.j
  }
  M.pr <- M/pr.iter
  
  ## Starting pooled pcm
  M <-0
  
  for (i in 1:r){ 
    ols1 <- dr.directions(dr(y[,i]~x, method="ols"))[,1]
    phd1 <- dr.directions(dr(y[,i]~x, method="phdres"))[,1]
    PI <- cbind(ols1, phd1)  
    
    if (nclust == 4) { kmeans.out <- dr.slices(PI, nslices=c(2,2))
    actual.k <- kmeans.out$nslices
    C.x <- kmeans.out$slice.indicator
    size.Cx <- kmeans.out$slice.sizes
    total.Cx <- sum(size.Cx) }  else { kmeans.out <- kmeans.min(PI, k=nclust, min=minc)
    actual.k <- kmeans.out$nslices
    C.x <- kmeans.out$slice.indicator
    size.Cx <- kmeans.out$slice.sizes
    total.Cx <- sum(size.Cx)   }
    
    
    ## Start coordinate mean
    M.j <- NULL
    
    for(k in 1:actual.k){
      sel.k <- C.x==k
      z.k <- z[sel.k, ] ;  y.k <- y[sel.k, i]
      wts.k <- wts[sel.k]
      
      slices <- dr.slices(y.k, nslices)
      zmeans <- matrix(0,slices$nslices,NCOL(z.k))
      slice.weight <- rep(0,slices$nslices)  # NOT rep(0,NCOL(z.k))
      
      # compute weighted means within slice 
      wmean <- function (x, wts) { sum(x * wts) / sum (wts) }
      for (j in 1:slices$nslices){
        sel <- slices$slice.indicator==j
        zmeans[j,]<- apply(z.k[sel, ],2, wmean, wts.k[sel])
        slice.weight[j]<-sum(wts.k[sel])}
      
      # get M matrix for sir
      ini.M <- t(apply(zmeans, 2, "*", sqrt(slice.weight))/ sqrt(sum(slice.weight)))
      M.j <- cbind(M.j, ( ini.M*(size.Cx[k]/total.Cx) ))
    } ## End coordinate mean
    
    M.j <- M.j%*% t(M.j)
    M <- M+M.j
  }## End for corrdinate regression
  
  M.pool <- M/r
  
  ## Summing PRpcm + pooled pcm
  M <- (M.pool+M.pr)*0.5
  
  return(list(M = M, slice.info = kmeans.out))
}

#######################################################################
## Multiveraite version of Reduced-Rank Response (RRR)
#######################################################################
dr.M.rrr<-function(object, nslices=2, dy=2, nclust=4, minc=5, slice.info=NULL,...) {
  z <- dr.z(object)
  x <- dr.x(object)
  y <- dr.y(object)
  wts <- dr.wts(object)
  
  r <- dim(y)[2]
  p <- dim(z)[2]
  
  ## Stating reduced-rank response mean method
  cand.mat <- NULL
  s0 <- mbrdr(y~x, method="prr")
  cand.mat <- cbind(cand.mat, s0$evectors[,1:dy])
  
  s1 <- update(s0, method="pfrr")
  cand.mat <- cbind(cand.mat, s1$evectors[,1:dy])
  
  s2 <- update(s0, method="upfrr")
  cand.mat <- cbind(cand.mat, s2$evectors[,1:dy])
  
  s3 <- update(s0, method="yc")
  cand.mat <- cbind(cand.mat, s3$evectors[,1:dy])
  
  pr.iter <- dim(cand.mat)[2]
  M<-0
  for (i in 1:pr.iter){ 
    t.i <- cand.mat[,i]
    t.i <- t.i /sqrt(sum(t.i^2))
    t.i.y <- c(y %*% t.i)
    ols1 <- dr.directions(dr( t.i.y~x, method="ols"))[,1]
    phd1 <- dr.directions(dr( t.i.y~x, method="phdres"))[,1]
    PI <- cbind(ols1, phd1)  
    
    if (nclust == 4) { kmeans.out <- dr.slices(PI, nslices=c(2,2))
    actual.k <- kmeans.out$nslices
    C.x <- kmeans.out$slice.indicator
    size.Cx <- kmeans.out$slice.sizes
    total.Cx <- sum(size.Cx) }  else { kmeans.out <- kmeans.min(PI, k=nclust, min=minc)
    actual.k <- kmeans.out$nslices
    C.x <- kmeans.out$slice.indicator
    size.Cx <- kmeans.out$slice.sizes
    total.Cx <- sum(size.Cx)   }
    
    M.j <- NULL
    for(k in 1:actual.k){
      sel.k <- C.x==k
      z.k <- z[sel.k, ] ;  y.k <- t.i.y[sel.k]
      wts.k <- wts[sel.k]
      
      slices <- dr.slices.arc(y.k, nslices)
      zmeans <- matrix(0,slices$nslices,NCOL(z.k))
      slice.weight <- rep(0,slices$nslices)  # NOT rep(0,NCOL(z.k))
      
      # compute weighted means within slice 
      wmean <- function (x, wts) { sum(x * wts) / sum (wts) }
      for (j in 1:slices$nslices){
        sel <- slices$slice.indicator==j
        zmeans[j,]<- apply(z.k[sel,],2, wmean, wts.k[sel])
        slice.weight[j]<-sum(wts.k[sel])}
      
      # get M matrix for sir
      ini.M <- t(apply(zmeans, 2, "*", sqrt(slice.weight))/ sqrt(sum(slice.weight)))
      M.j <- cbind(M.j, ( ini.M*(size.Cx[k]/total.Cx) ))
    } ## coordinate mean
    
    ##    if (dim(M.j)[2] != (nslices*nclust) ) {
    ##         diff.col <-  (nslices*nclust) - dim(M.j)[2]
    ##         M.j.zeros <- matrix( rep(0, (p*diff.col) ), c(p, diff.col) )
    ##         M.j <- cbind(M.j, M.j.zeros) }
    
    M.j <- M.j%*% t(M.j)
    M <- M+M.j
  }
  M.pr <- M/pr.iter
  
  return(list(M = M.pr, slice.info = kmeans.out))
}

##
## Multiveraite version of PCM
##
#######################################################################
#######################################################################

#######################################################################
## Pooled PCM + Projection Resampling PCM
#######################################################################
dr.M.PRpoolpcm<-function(object, nslices=2, nclust=4, pr.iter=100, minc=5, slice.info=NULL,...) {
  z <- dr.z(object)
  x <- dr.x(object)
  y <- dr.y(object)
  wts <- dr.wts(object)
  
  r <- dim(y)[2]
  p <- dim(z)[2]
  
  ## Stating PRpcm
  M <-0
  
  for (i in 1:pr.iter){ 
    t.i <- as.matrix(rnorm(r))
    t.i <- t.i /sqrt(sum(t.i^2))
    t.i.y <- c(y %*% t.i)
    ols1 <- dr.directions(dr( t.i.y~x, method="ols"))[,1]
    phd1 <- dr.directions(dr( t.i.y~x, method="phdres"))[,1]
    PI <- cbind(ols1, phd1)  
    
    if (nclust == 4) { kmeans.out <- dr.slices(PI, nslices=c(2,2))
    actual.k <- kmeans.out$nslices
    C.x <- kmeans.out$slice.indicator
    size.Cx <- kmeans.out$slice.sizes
    total.Cx <- sum(size.Cx) }  else { kmeans.out <- kmeans.min(PI, k=nclust, min=minc)
    actual.k <- kmeans.out$nslices
    C.x <- kmeans.out$slice.indicator
    size.Cx <- kmeans.out$slice.sizes
    total.Cx <- sum(size.Cx)   }
    
    M.j <- NULL
    for(k in 1:actual.k){
      sel.k <- C.x==k
      z.k <- z[sel.k, ] ;  y.k <- t.i.y[sel.k]
      wts.k <- wts[sel.k]
      
      slices <- dr.slices.arc(y.k, nslices)
      zmeans <- matrix(0,slices$nslices,NCOL(z.k))
      slice.weight <- rep(0,slices$nslices)  # NOT rep(0,NCOL(z.k))
      
      # compute weighted means within slice 
      wmean <- function (x, wts) { sum(x * wts) / sum (wts) }
      for (j in 1:slices$nslices){
        sel <- slices$slice.indicator==j
        zmeans[j,]<- apply(z.k[sel,],2, wmean, wts.k[sel])
        slice.weight[j]<-sum(wts.k[sel])}
      
      # get M matrix for sir
      ini.M <- t(apply(zmeans, 2, "*", sqrt(slice.weight))/ sqrt(sum(slice.weight)))
      M.j <- cbind(M.j, ( ini.M*(size.Cx[k]/total.Cx) ))
    } ## coordinate mean
    
    ##    if (dim(M.j)[2] != (nslices*nclust) ) {
    ##         diff.col <-  (nslices*nclust) - dim(M.j)[2]
    ##         M.j.zeros <- matrix( rep(0, (p*diff.col) ), c(p, diff.col) )
    ##         M.j <- cbind(M.j, M.j.zeros) }
    
    M.j <- M.j%*% t(M.j)
    M <- M+M.j
  }
  M.pr <- M/pr.iter
  
  ## Starting pooled pcm
  M <-0
  
  for (i in 1:r){ 
    ols1 <- dr.directions(dr(y[,i]~x, method="ols"))[,1]
    phd1 <- dr.directions(dr(y[,i]~x, method="phdres"))[,1]
    PI <- cbind(ols1, phd1)  
    
    if (nclust == 4) { kmeans.out <- dr.slices(PI, nslices=c(2,2))
    actual.k <- kmeans.out$nslices
    C.x <- kmeans.out$slice.indicator
    size.Cx <- kmeans.out$slice.sizes
    total.Cx <- sum(size.Cx) }  else { kmeans.out <- kmeans.min(PI, k=nclust, min=minc)
    actual.k <- kmeans.out$nslices
    C.x <- kmeans.out$slice.indicator
    size.Cx <- kmeans.out$slice.sizes
    total.Cx <- sum(size.Cx)   }
    
    
    ## Start coordinate mean
    M.j <- NULL
    
    for(k in 1:actual.k){
      sel.k <- C.x==k
      z.k <- z[sel.k, ] ;  y.k <- y[sel.k, i]
      wts.k <- wts[sel.k]
      
      slices <- dr.slices(y.k, nslices)
      zmeans <- matrix(0,slices$nslices,NCOL(z.k))
      slice.weight <- rep(0,slices$nslices)  # NOT rep(0,NCOL(z.k))
      
      # compute weighted means within slice 
      wmean <- function (x, wts) { sum(x * wts) / sum (wts) }
      for (j in 1:slices$nslices){
        sel <- slices$slice.indicator==j
        zmeans[j,]<- apply(z.k[sel, ],2, wmean, wts.k[sel])
        slice.weight[j]<-sum(wts.k[sel])}
      
      # get M matrix for sir
      ini.M <- t(apply(zmeans, 2, "*", sqrt(slice.weight))/ sqrt(sum(slice.weight)))
      M.j <- cbind(M.j, ( ini.M*(size.Cx[k]/total.Cx) ))
    } ## End coordinate mean
    
    M.j <- M.j%*% t(M.j)
    M <- M+M.j
  }## End for corrdinate regression
  
  M.pool <- M/r
  
  ## Summing PRpcm + pooled pcm
  M <- (M.pool+M.pr)*0.5
  
  return(list(M = M, slice.info = kmeans.out))
}



#######################################################################
## Pooled PCM
#######################################################################
dr.M.poolpcm<-function(object, nslices=2, nclust=4, minc=10, slice.info=NULL,...) {
  z <- dr.z(object);   x <- dr.x(object)
  y <- dr.y(object);    wts <- dr.wts(object)
  
  r <- dim(y)[2];     p <- dim(z)[2]
  
  M <-0
  
  for (i in 1:r){ 
    ols1 <- dr.directions(dr(y[,i]~x, method="ols"))[,1]
    phd1 <- dr.directions(dr(y[,i]~x, method="phdres"))[,1]
    PI <- cbind(ols1, phd1)  
    
    if (nclust == 4) { kmeans.out <- dr.slices(PI, nslices=c(2,2))
    actual.k <- kmeans.out$nslices
    C.x <- kmeans.out$slice.indicator
    size.Cx <- kmeans.out$slice.sizes
    total.Cx <- sum(size.Cx) }  else { kmeans.out <- kmeans.min(PI, k=nclust, min=minc)
    actual.k <- kmeans.out$nslices
    C.x <- kmeans.out$slice.indicator
    size.Cx <- kmeans.out$slice.sizes
    total.Cx <- sum(size.Cx)   }
    
    
    ## Start coordinate mean
    M.j <- NULL
    
    for(k in 1:actual.k){
      sel.k <- C.x==k
      z.k <- z[sel.k, ] ;  y.k <- y[sel.k, i]
      wts.k <- wts[sel.k]
      
      slices <- dr.slices(y.k, nslices)
      zmeans <- matrix(0,slices$nslices,NCOL(z.k))
      slice.weight <- rep(0,slices$nslices)  # NOT rep(0,NCOL(z.k))
      
      # compute weighted means within slice 
      wmean <- function (x, wts) { sum(x * wts) / sum (wts) }
      for (j in 1:slices$nslices){
        sel <- slices$slice.indicator==j
        zmeans[j,]<- apply(z.k[sel, ],2, wmean, wts.k[sel])
        slice.weight[j]<-sum(wts.k[sel])}
      
      # get M matrix for sir
      ini.M <- t(apply(zmeans, 2, "*", sqrt(slice.weight))/ sqrt(sum(slice.weight)))
      M.j <- cbind(M.j, ( ini.M*(size.Cx[k]/total.Cx) ))
    } ## End coordinate mean
    
    M.j <- M.j%*% t(M.j)
    M <- M+M.j
  }## End for corrdinate regression
  
  M <- M / r
  
  return(list(M = M, slice.info = kmeans.out))
}

#######################################################################
## Projection Resampling PSM
#######################################################################
dr.M.PRpcm<-function(object, nslices=2, nclust=4, pr.iter=1000, minc=5, slice.info=NULL,...) {
  z <- dr.z(object)
  x <- dr.x(object)
  y <- dr.y(object)
  wts <- dr.wts(object)
  
  r <- dim(y)[2]
  p <- dim(z)[2]
  
  M <-0
  
  for (i in 1:pr.iter){ 
    t.i <- as.matrix(rnorm(r))
    t.i <- t.i /sqrt(sum(t.i^2))
    t.i.y <- c(y %*% t.i)
    ols1 <- dr.directions(dr( t.i.y~x, method="ols"))[,1]
    phd1 <- dr.directions(dr( t.i.y~x, method="phdres"))[,1]
    PI <- cbind(ols1, phd1)  
    
    if (nclust == 4) { kmeans.out <- dr.slices(PI, nslices=c(2,2))
    actual.k <- kmeans.out$nslices
    C.x <- kmeans.out$slice.indicator
    size.Cx <- kmeans.out$slice.sizes
    total.Cx <- sum(size.Cx) }  else { kmeans.out <- kmeans.min(PI, k=nclust, min=minc)
    actual.k <- kmeans.out$nslices
    C.x <- kmeans.out$slice.indicator
    size.Cx <- kmeans.out$slice.sizes
    total.Cx <- sum(size.Cx)   }
    
    M.j <- NULL
    for(k in 1:actual.k){
      sel.k <- C.x==k
      z.k <- z[sel.k, ] ;  y.k <- t.i.y[sel.k]
      wts.k <- wts[sel.k]
      
      slices <- dr.slices.arc(y.k, nslices)
      zmeans <- matrix(0,slices$nslices,NCOL(z.k))
      slice.weight <- rep(0,slices$nslices)  # NOT rep(0,NCOL(z.k))
      
      # compute weighted means within slice 
      wmean <- function (x, wts) { sum(x * wts) / sum (wts) }
      for (j in 1:slices$nslices){
        sel <- slices$slice.indicator==j
        zmeans[j,]<- apply(z.k[sel,],2, wmean, wts.k[sel])
        slice.weight[j]<-sum(wts.k[sel])}
      
      # get M matrix for sir
      ini.M <- t(apply(zmeans, 2, "*", sqrt(slice.weight))/ sqrt(sum(slice.weight)))
      M.j <- cbind(M.j, ( ini.M*(size.Cx[k]/total.Cx) ))
    } ## coordinate mean
    
    ##    if (dim(M.j)[2] != (nslices*nclust) ) {
    ##         diff.col <-  (nslices*nclust) - dim(M.j)[2]
    ##         M.j.zeros <- matrix( rep(0, (p*diff.col) ), c(p, diff.col) )
    ##         M.j <- cbind(M.j, M.j.zeros) }
    
    M.j <- M.j%*% t(M.j)
    M <- M+M.j
  }
  
  M <- M/pr.iter
  
  return(list(M = M, slice.info = kmeans.out))
}

# Model 1

n<-100
r42 <-NULL
for (i in 1:500){
  u1 <- runif(n); e <- runif(n,-0.5,0.5); u2 <- log(u1)+e
  w1 <- rnorm(n); w2 <- rnorm(n); w3 <- rnorm(n)
  
  x1 <- u1+w1 ; x2 <- u2+w2+w3; x3 <- w1-w2; x4 <- w2; x5 <- w3
  x6 <- rnorm(n); x7 <- rnorm(n); x8 <- rnorm(n)
  x9 <- rnorm(n) ; x10 <- rnorm(n)
  
  eta.x <- x1-x2-x3
  y1 <- exp(0.5*eta.x+1) +0.1*rnorm(n)
  y2 <- eta.x^2 +0.1*rnorm(n)
  y3 <- y1+y2 +0.1*rnorm(n)
  y4 <-abs(eta.x) +0.1*rnorm(n)
  
  sim.dat<-data.frame(y1, y2, y3, y4, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)
  
  rrr <- dr(cbind(y1, y2, y3, y4)~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data=sim.dat, dy=2, nclust=4, nslices=2, method="rrr")
  
  eta.hat.x <- dr.directions(rrr)[,1]
  
  r42[i]<- sqrt( summary(lm(eta.x~eta.hat.x))$r.squared)
} 



