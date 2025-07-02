get_Init<-function(M,U,u,G=NULL){ 
  MU = M+U
  invMU<- chol2inv(chol(MU))
  obj<-function(G){
    logdet(crossprod(G,M)%*%G)+logdet(crossprod(G,invMU)%*%G)
  }
  
  # get init1  
  tmp.MU <- eigen(MU)
  tmp2.MU <- crossprod(tmp.MU$vectors, U) %*% tmp.MU$vectors
  tmp3.MU <- sort(diag(tmp2.MU), decreasing = TRUE, index.return = TRUE)
  init <- as.matrix(tmp.MU$vectors[, tmp3.MU$ix[1:u]])   
  obj1 <- obj(init)
  if(!missing(G)) print(subspace(init,G))
  
  # get init2	
  tmp2 <- diag(1/sqrt(tmp.MU$values)) 
  tmp3 <- sort(diag(tmp2 %*% tmp2.MU %*% tmp2), decreasing = TRUE, index.return = TRUE)
  init.MU <- as.matrix(tmp.MU$vectors[, tmp3$ix[1:u]])
  if(!missing(G))print(subspace(init.MU,G))
  
  obj2 <- obj(init.MU)
  if (obj2 < obj1) {
    init <- init.MU
    obj1 <- obj2
  }
  
  # get init3
  tmp.M <- eigen(M)
  tmp2.M <- crossprod(tmp.M$vectors, U) %*% tmp.M$vectors
  tmp3.M <- sort(diag(tmp2.M), decreasing = TRUE, index.return = TRUE)	
  init.M <- as.matrix(tmp.M$vectors[, tmp3.M$ix[1:u]])
  obj3 <-  obj(init.M)
  if (obj3 < obj1) {
    init <- init.M
    obj1 <- obj3
  }
  if(!missing(G)) print(subspace(init.M,G))
  
  # get init4				
  tmp2.M <- diag(1/sqrt(tmp.M$values)) 
  tmp3.M <- sort(diag(tmp2 %*% tmp2.M %*% tmp2), decreasing = TRUE, index.return = TRUE)
  init.M <- as.matrix(tmp.M$vectors[, tmp3.M$ix[1:u]])
  if(!missing(G)) print(subspace(init.M,G))
  obj4 <-  obj(init.M)
  if (obj4 < obj1) {
    init <- init.M
    obj1 <- obj4
  }      
  return(init) 
}


