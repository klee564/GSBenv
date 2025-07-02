initial_value3 <- function(X, Y, u,verbose=0,choose=NULL) {
  X <-as.matrix(X)
  Y <-as.matrix(Y)
  n <- nrow(Y)
  r <- ncol(Y)
  p <- ncol(X)
  #mX = colMeans(X)
  #mY = colMeans(Y)
  Yc <- as.matrix(scale(Y, center = T, scale = FALSE))
  Xc <- as.matrix(scale(X, center = T, scale = FALSE))
  sigY <- cov(Yc)*(n-1)/n
  sigX <- cov(Xc)*(n-1)/n 
  sigYX <- cov(Yc, Xc)*(n-1)/n
  if(n>r) invsigY <- chol2inv(chol(sigY))
  
  tmp = chol(sigX) # t(tmp)%*%tmp =sigX
  invsigX <- chol2inv(tmp) # invsigX=tmp2%*%t(tmp2)
  tmp2 <- backsolve(tmp, diag(p)) # tmp2=inv(tmp)  
  tmp3 <- sigYX %*% tmp2
  bsxb <- tcrossprod(tmp3,tmp3) #sigYX %*% invsigX%*% t(sigYX)
  sigRes <- sigY - bsxb
  
  obj<-function(G){
    tmp1=crossprod(G,sigRes)%*%G
    tmp2=crossprod(G,invsigY)%*%G
    logdet(tmp1)+logdet(tmp2)
  }
  
  #=============
  
  tmp.y <- eigen(sigY)
  tmp2.y <- crossprod(tmp.y$vectors, bsxb) %*% tmp.y$vectors
  tmp3.y <- sort(diag(tmp2.y), decreasing = TRUE, index.return = TRUE)
  ##print(sort(tmp3.y$ix[1:u]))
  init <- as.matrix(tmp.y$vectors[, tmp3.y$ix[1:u]]) 
  #if(verbose) {cat('init1\n');print(init)}
  k <- 1	
  #print(init)	
  if (n > r + 1) {
    obj1 <- obj(init)			
    tmp2 <- diag(1 / sqrt(tmp.y$values))
    tmp3 <- sort(diag(tmp2 %*% tmp2.y %*% tmp2), decreasing = TRUE, index.return = TRUE)
    ##print(sort(tmp3$ix[1:u]))
    init.y <- as.matrix(tmp.y$vectors[, tmp3$ix[1:u]])
    #if(verbose) {cat('init2\n');print(init.y)}
    obj2 <- sum(obj(init.y))		
    if (obj2 < obj1) {
      init <- init.y
      obj1 <- obj2
      k <- 2
    }
    if (n > r + p + 1) {
      tmp.res <- eigen(sigRes)
      tmp2.res <- crossprod(tmp.res$vectors, bsxb) %*% tmp.res$vectors
      tmp3.res <- sort(diag(tmp2.res), decreasing = TRUE, index.return = TRUE)	
      init.res <- as.matrix(tmp.res$vectors[, tmp3.res$ix[1:u]])				
      #if(verbose) {cat('init3\n');print(init.res)}
      obj3 <- obj(init.res)		
      if (obj3 < obj1) {
        init <- init.res
        obj1 <- obj3
        k <- 3
      }
      
      tmp2 <- diag(1 / sqrt(tmp.res$values)) 
      tmp3.res <- sort(diag(tmp2 %*% tmp2.res %*% tmp2), decreasing = TRUE, index.return = TRUE)
      init.res2 <- as.matrix(tmp.res$vectors[, tmp3.res$ix[1:u]])
      #if(verbose) {cat('init4\n');print(init.res)}
      
      obj4 <- obj(init.res2)
      if (obj4 < obj1) {
        init <- init.res2
        obj1 <- obj4
        k <- 4
      }
    }
  }
  if(verbose)	cat('We choose init',k,'\n')		
  result=list(init1=init,init2=init.y,init3= init.res, init4=  init.res2, objs=c(obj1,obj2,obj3,obj4))
  return(result)
}



