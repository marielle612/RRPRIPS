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

