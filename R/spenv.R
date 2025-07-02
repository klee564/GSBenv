spenv <- function(X, Y, u, lambda1=NULL,lambda2=NULL, ftol = 1e-4, maxiter = 1e2,verbose=0,init=NULL) {
  t1 = proc.time()
  X = as.matrix(X)
  Y = as.matrix(Y)
  r = ncol(Y)
  if(u==0|u==r) {out=env(X,Y,u)}
  else{
    if(is.null(lambda1)) lambda1 <- exp(seq(log(1),log(1e-5),len = 15))
    if(is.null(lambda2)) lambda2 <- exp(seq(log(1),log(1e-5),len = 15))
    if(is.null(init))  init =  initial_value(X,Y,u)
    GEidx = GE(init)
    newY = Y[, GEidx]
    newinit = init[GEidx,,drop = FALSE]

    m1 <- LassoLambda.spenv(X, newY, u, lambda=lambda1, ftol = ftol,
                            maxiter = maxiter, weight = rep(1,r-u),
                            init = newinit,verbose=verbose)
    #calculating weight
    Gammahat <- m1$Gamma
    w <- Gammahat %*% solve(Gammahat[1:u, ])
    w_norm <- 1/(rowSums(w^2)^2)[(u+1):r]
    #w_norm2 <- 1/(rowSums(Gammahat^2)^2)[(u+1):r]
    #print(w_norm/w_norm2)

    if(verbose) print(w_norm)
    #adaptive lasso step
    #cat('-----------Second Adaptive Lasso\n')
    #newinit = m1$Gamma;
    m2 <- LassoLambda.spenv(X, newY, u, lambda=lambda2, ftol = ftol,
                            maxiter = maxiter, weight = w_norm,
                            init = newinit,verbose=verbose)
    out = m2
    out$alpha = m2$alpha[order(GEidx),,drop=FALSE]
    out$Gamma = m2$Gamma[order(GEidx),,drop=FALSE]
    out$Gamma0 = m2$Gamma0[order(GEidx),,drop=FALSE]
    out$beta = m2$beta[order(GEidx),,drop=FALSE]
    out$where1 =  sort(GEidx[m2$where1])
    out$where0 =  setdiff(1:r,out$where1)
    out$BIC_seq = NULL
    out$BIC_seq1=m1$BIC_seq
    out$BIC_seq2=m2$BIC_seq
    out$BIC=m2$BIC
    out$sigRes=m2$sigRes[order(GEidx),order(GEidx)]
    out$sigY=m2$sigY[order(GEidx),order(GEidx)]
    out$lambda=c(m1$lambda,m2$lambda)
  }
  out$fit_time = (proc.time()-t1)[3]
  class(out)<-'spenv'
  out
}



summary_senv <- function(senv.out){
  r <- nrow(senv.out$beta)
  xi <- rep(1,r)
  xi[senv.out$where0] <-0

  list(beta= senv.out$beta, xi=xi)
}
