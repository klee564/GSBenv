get_Init2<-function(M,U,u){ 
  MU = M+U
  #invMU<- chol2inv(chol(MU))
  tmp.MU <- eigen(MU)
  tmp2 <- diag(1/sqrt(tmp.MU$values)) 
  tmp.M <- eigen(M)

  tmp2.M <- diag(1/sqrt(tmp.M$values)) 
  tmp3.M <- sort(diag(tmp2 %*% tmp2.M %*% tmp2), decreasing = TRUE, index.return = TRUE)
  init.M <- as.matrix(tmp.M$vectors[, tmp3.M$ix[1:u]])
  
  
  return(init.M) 
}


