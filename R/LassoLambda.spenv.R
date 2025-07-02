LassoLambda.spenv <- function(X, Y, u, lambda=NULL, ftol=1e-2, maxiter=1e2, weight = NULL,init=NULL,verbose=0) {
  t1 = proc.time()
  X = as.matrix(X)
  Y = as.matrix(Y)
  r = ncol(Y)
  n = nrow(Y)
  if(missing(weight)) weight=rep(1,r-u)
  if(missing(lambda)) lambda <- exp(seq(log(1),log(1e-5),len=15))


  if(missing(init)) init <- initial_value(X,Y,u)
  else init=init

  model_vec <- rep(NA, length(lambda))
  names(model_vec)=signif(lambda,2)
  minBIC = Inf
  res=NULL
  idmin = NULL
  for(l in 1:length(lambda)){
    if(verbose) cat('----------------lambda =',lambda[l],'\n')
    tmp <- spenvbase(X, Y, u,	ftol=ftol, 	maxiter=maxiter, lambda=lambda[l], weight=weight,init=init,verbose=verbose)
  	init=tmp$Gamma #warm start
    BIC <- -2*tmp$loglik + log(tmp$n) * (tmp$q-u) * u
    model_vec[l] <- BIC
    if(BIC<minBIC) {minBIC=BIC
                    res=tmp
                    idmin=l
    }
  }

  lambda.min <- lambda[idmin]
  out <- res
  out$lambda=lambda.min
  out$BIC_seq=model_vec
  out$fit_time=(proc.time()-t1)[3]
  out$BIC <- model_vec[idmin]
  out
}
