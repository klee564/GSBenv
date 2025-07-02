spenvbase <- function(X, Y, u,lambda, weight=NULL, ftol=1e-4, eps2=1e-4,maxiter=1e2,init=NULL,verbose=0) {
  t1 = proc.time()
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  n <- nrow(Y)
  r <- ncol(Y)
  p <- ncol(X)
  if(missing(weight)) weight=rep(1,r-u)
  mX = colMeans(X)
  mY = colMeans(Y)
  Yc <- as.matrix(scale(Y, center = T, scale = FALSE))
  Xc <- as.matrix(scale(X, center = T, scale = FALSE))
  sigY <- cov(Yc)*(n-1)/n
  sigX <- cov(Xc)*(n-1)/n
  sigYX <- cov(Yc, Xc)*(n-1)/n   
  invsigY <- chol2inv(sechol(sigY))   
  tmp = sechol(sigX) # t(tmp)%*%tmp =sigX
  invsigX <- chol2inv(tmp) # invsigX=tmp2%*%t(tmp2)
  tmp2 <- backsolve(tmp, diag(p)) # tmp2=inv(tmp)  
  tmp3 <- sigYX %*% tmp2
  U <- tcrossprod(tmp3,tmp3) #sigYX %*% invsigX%*% t(sigYX)
  sigRes <- sigY - U
  betaOLS <- sigYX %*% invsigX  
  logDetsigY = logdet(sigY)
  ModelOutput=list()
  
  if (u == 0) {    
    ModelOutput$alpha = mY;
    ModelOutput$beta = matrix(0,r, p)
    ModelOutput$Gamma = NULL
    ModelOutput$Gamma0 = diag(r);
    ModelOutput$eta = NULL;
    ModelOutput$Sigma = sigY    
    ModelOutput$Omega = NULL;
    ModelOutput$Omega0 = sigY;
    ModelOutput$loglik = - n * r / 2 * (1 + log(2 * pi)) - n / 2 * logDetsigY;
    ModelOutput$paramNum = r + u * p + r * (r + 1) / 2;
    ModelOutput$sigRes=sigRes   
    ModelOutput$sigY=sigY
    ModelOutput$sigX=sigX
    ModelOutput$n = n 
    ModelOutput$r = r
    ModelOutput$u = u
    ModelOutput$p = p
    ModelOutput$q = 0
    ModelOutput$where1 = NULL
    ModelOutput$where0 = 1:r
    ModelOutput$lambda = lambda
  } 
  else if (u == r) {    
    ModelOutput$alpha = mY - betaOLS * mX;
    ModelOutput$beta = betaOLS;
    ModelOutput$Gamma = diag(r);
    ModelOutput$Gamma0 = NULL;
    ModelOutput$eta = betaOLS;
    ModelOutput$Sigma = sigRes;
    ModelOutput$Omega = sigRes;
    ModelOutput$Omega0 = NULL;
    ModelOutput$loglik = - n * r / 2 * (1 + log(2 * pi)) - n / 2 * logdet(sigRes);
    ModelOutput$paramNum = r + u * p + r * (r + 1) / 2;
    ModelOutput$sigRes=sigRes   
    ModelOutput$sigY=sigY
    ModelOutput$sigX=sigX
    ModelOutput$n = n
    ModelOutput$r = r
    ModelOutput$u = u
    ModelOutput$p = p
    ModelOutput$q = r
    ModelOutput$where1 = 1:r
    ModelOutput$where0 = NULL
    ModelOutput$lambda=lambda
  } 
  else{
    if(missing(init)) init <- initial_value(X,Y,u)
    Ginit <- init %*% solve(init[1:u, ])
    if (u == (r-1)) {
      obj1 <- logdet(t(init) %*% sigRes %*% init) + logdet(t(init) %*% invsigY %*% init)
      widx <- is.finite(weight)      
      Ginit.tmp <- as.matrix(Ginit[(u+1):r, ])
      Ginit.tmp <- as.matrix(Ginit.tmp[widx,])
      obj1 <- obj1 + lambda * sum(weight[widx] * sqrt(rowSums(Ginit.tmp^2)))
      obj_seq = obj1
      
      U1c2 <- array(0, dim = c(r-1, r-1))
      V1c2 <- array(0, dim = c(r-1, r-1))      
      U1c2 <- sigRes[-r, -r] - as.matrix(sigRes[-r, r]) %*% sigRes[r, -r] / sigRes[r, r]
      V1c2 <- invsigY[-r, -r] - as.matrix(invsigY[-r, r]) %*% invsigY[r, -r] / invsigY[r, r]	      
      t2 <- sigRes[-r, r] / sigRes[r, r]
      t3 <- invsigY[-r, r] / invsigY[r, r]
      invC1 <- chol2inv(sechol(U1c2))
      invC2 <- chol2inv(sechol(V1c2))
      i <- 1		
      while (i < maxiter) {
        res <- spenvlp(b2=drop(t2), b3=drop(t3),A1=diag(r-1), 
                       A2=sigRes[r, r]*invC1, A3=invsigY[r, r]*invC2,
                       lambda=lambda,eps=eps2, maxiter=maxiter,
                       weight=weight[r-u],a_vec_init=drop(Ginit[r,]))      
        #old_Ginit <- Ginit[r, ]
        Ginit[r, ] <- res$a_vec
        
        a <- qr(Ginit)
        Gamma <- qr.Q(a)
        Ginit.tmp <- as.matrix(Ginit[(u+1):r, ])
        Ginit.tmp <- as.matrix(Ginit.tmp[widx,])
        obj5 <- logdet(t(Gamma) %*% sigRes %*% Gamma)+ logdet(t(Gamma) %*% invsigY %*% Gamma) + lambda * sum(weight[widx] * sqrt(rowSums(Ginit.tmp^2)))
        obj_seq <- c(obj_seq , obj5)
        #print('------------')
        if (abs(obj1 - obj5) < ftol * abs(obj1)) {
          if(verbose) print(obj_seq)
          break
        }	
        else {
          obj1 <- obj5
          i <- i + 1
        }
        #if(sum((Ginit[r,]-old_Ginit)^2) < 1e-4) break
        #i <- i + 1		
      }
      Gamma=qr.Q(qr(Ginit))
    } 
    else {        
      obj1 <- logdet(t(init) %*% sigRes %*% init) + logdet(t(init) %*% invsigY %*% init)
      widx <- is.finite(weight)      
      Ginit.tmp <- as.matrix(Ginit[(u+1):r, ])
      Ginit.tmp <- as.matrix(Ginit.tmp[widx,])
      obj1 <- obj1 + lambda * sum(weight[widx] * sqrt(rowSums(Ginit.tmp^2)))
      obj_seq = obj1
      
      GUG <- crossprod(Ginit, (sigRes %*% Ginit))	
      GVG <- crossprod(Ginit, (invsigY %*% Ginit))		
      
      t4 <- crossprod(Ginit[(u+1):r,], Ginit[(u+1):r, ]) + diag(u)
      i <- 1
      while (i < maxiter) {
        #old_Ginit <- Ginit
        for (j in (u+1):r) {
          g <- as.matrix(Ginit[j, ])
          t2 <- crossprod(Ginit[-j, ], as.matrix(sigRes[-j, j])) / sigRes[j, j]
          t3 <- crossprod(Ginit[-j, ], as.matrix(invsigY[-j, j])) / invsigY[j, j]
          
          GUGt2 <- g + t2
          GUG <- GUG - tcrossprod(GUGt2, GUGt2) * sigRes[j, j]
          
          GVGt2 <- g + t3
          GVG <- GVG - tcrossprod(GVGt2, GVGt2) * invsigY[j, j] 
          
          t4 <- t4 - tcrossprod(as.matrix(Ginit[j, ]), as.matrix(Ginit[j, ]))
          invC1 <- chol2inv(sechol(GUG))
          invC2 <- chol2inv(sechol(GVG))
          invt4 <- chol2inv(sechol(t4))
          
          res <- spenvlp(b2=drop(t2), b3=drop(t3), A1=invt4, A2=sigRes[j, j]*invC1, A3=invsigY[j, j]*invC2, lambda=lambda, eps=eps2, maxiter=maxiter, weight=weight[j-u], a_vec_init=drop(Ginit[j,]))
          #           if(i==1&&j==u+1){
          #             #print(res$vecs)
          #             innerloop_objseq=rep(0,nrow(res$vecs))
          #             for(l in 1:nrow(res$vecs)){
          #               Ginit[j, ]=res$vecs[l,]
          #               a <- qr(Ginit)
          #               Gamma <- qr.Q(a)
          #               Ginit.tmp <- as.matrix(Ginit[(u+1):r, ])
          #               Ginit.tmp <- as.matrix(Ginit.tmp[widx,])
          #               innerloop_objseq[l] <- logdet(t(Gamma) %*% sigRes %*% Gamma)+ logdet(t(Gamma) %*% invsigY %*% Gamma) + lambda * sum(weight[widx] * sqrt(rowSums(Ginit.tmp^2)))
          #             }
          #             print(innerloop_objseq)
          #}
          Ginit[j, ] <- res$a_vec
          g <- as.matrix(Ginit[j, ])
          t4 <- t4 + tcrossprod(g, g)
          GUGt2 <- g + t2
          GUG <- GUG + tcrossprod(GUGt2, GUGt2) * sigRes[j, j]
          
          GVGt2 <- g + t3
          GVG <- GVG + tcrossprod(GVGt2, GVGt2) * invsigY[j, j] 
        }
        a <- qr(Ginit)
        Gamma <- qr.Q(a)
        Ginit.tmp <- as.matrix(Ginit[(u+1):r, ])
        Ginit.tmp <- as.matrix(Ginit.tmp[widx,])
        obj5 <- logdet(t(Gamma) %*% sigRes %*% Gamma)+ logdet(t(Gamma) %*% invsigY %*% Gamma) + lambda * sum(weight[widx] * sqrt(rowSums(Ginit.tmp^2)))
        obj_seq <- c(obj_seq , obj5)
        #print('------------')
        if (abs(obj1 - obj5) < ftol * abs(obj1)) {
          #if(verbose) print(obj_seq)
          break
        }	
        else {
          obj1 <- obj5
          i <- i + 1
        }		
      }
      if(verbose) cat('The number of iterations:',i,'.\n',sep='')     
    }
    Gamma0 = grams(nulbasis(t(Gamma)))
    eta <- crossprod(Gamma, betaOLS)
    beta <- Gamma %*% eta
    alpha = mY - beta %*% mX  
    Omega = t(Gamma) %*% sigRes %*% Gamma;
    Omega0 = t(Gamma0) %*% sigY %*% Gamma0
    Sigma1 = Gamma %*% Omega %*% t(Gamma)
    Sigma2 = Gamma0 %*% Omega0 %*% t(Gamma0)
    Sigma = Sigma1 + Sigma2   
    if(!is.na(sum(Gamma))) idx <- which((rowSums(abs(Gamma))>0))
    q=length(idx)
    paramNum = r + u * (p - r + q) + r * (r + 1) / 2  
    a = logdet(t(Gamma) %*% sigRes %*% Gamma)
    b = logdet(t(Gamma) %*% invsigY %*% Gamma)
    l = n * r * (1 + log(2 * pi)) + n * (a + b + logDetsigY)
    
    ModelOutput$alpha = alpha;
    ModelOutput$beta = beta;
    ModelOutput$Gamma = Gamma;
    ModelOutput$Gamma0 = Gamma0;
    ModelOutput$eta = eta;
    ModelOutput$Sigma = Sigma
    ModelOutput$Omega = Omega
    ModelOutput$Omega0 = Omega0
    ModelOutput$loglik = -0.5* l
    ModelOutput$paramNum = paramNum;
    ModelOutput$sigRes=sigRes   
    ModelOutput$sigY=sigY
    ModelOutput$sigX=sigX
    ModelOutput$obj_seq=obj_seq
    ModelOutput$n = n
    ModelOutput$r = r
    ModelOutput$u = u
    ModelOutput$p = p
    ModelOutput$q = q
    ModelOutput$where1 = idx
    ModelOutput$where0 = setdiff(1:r,idx)
    ModelOutput$lambda=lambda
    #ModelOutput$innerloop_objseq=innerloop_objseq
  }# end of else 
  ModelOutput$fit_time=(proc.time()-t1)[3]
  class(ModelOutput)<-'spenvbase'
  ModelOutput 
}






