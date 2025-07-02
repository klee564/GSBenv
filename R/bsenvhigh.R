
#'
#' @examples
#' r <- 5
#' p <- 4
#' u <- 2
#' xi.tru <- c(1,1,0)
#'
#' groupind <- c(1,1,2)
#' K <- max(groupindex)
#' Betashapes <- matrix(1/1000,K,2)
#' group_prior <- list(groupind=groupind,Betashapes=Betashapes)
#'
#' param <- generate_sparsepar(r,p,u, xi.tru)
#'
#' inputdata <- generate_data(param,200)
#' X=inputdata$X;Y=inputdata$Y
#' postlist <- envGS_MH(X=inputdata$X,Y=inputdata$Y,u,group_prior,tau0=0.01,tau1=1)
#' res <- purrr::map(postlist,~param_invtrans(.x$mu,.x$etatilde,.x$Omegatilde,.x$Omega0tilde,.x$A))
#' mean(((res %>% purrr::map(~.x$beta) %>% purrr::reduce(`+`))/100- param$beta.tru)^2)
#' mean((Renvlp::env(X=inputdata$X,Y=inputdata$Y,u = u)$beta- param$beta.tru)^2)
#'
#' mean(((res %>% purrr::map(~.x$mu) %>% purrr::reduce(`+`))/100- param$mu.tru)^2)
#' mean((Renvlp::env(X=inputdata$X,Y=inputdata$Y,u = u)$mu- param$mu.tru)^2)
#'
envGS_MH_high <- function(X,Y,u,group_prior=group_prior,tau0,tau1,mcmc.num=20,trueparam){

  postlist <- list()

  r <- ncol(Y)
  p <- ncol(X)
  n <- nrow(X)

  groupindex <- group_prior$groupind
  K <- max(groupindex)
  phi <- rep(1/2,K)

  barX <- colMeans(X)
  barY <- colMeans(Y)

  #set hyperparameter
  M <- diag(10^(-6),p)
  Msqrt <- sqrtmat(M)
  B0 <- matrix(0,r,p)
  psi <- 1
  psi0 <- 1
  nu <- u +1
  nu0 <- r-u +1

  Betashapes <- group_prior$Betashapes


  #set init
  init <- env.high(X=X,Y=Y,u = u)
  initA <- as.matrix(find_A_from_gamma(init$Gamma),col=u)
  gamma_gamma0 <- find_gammas_from_A(initA)
  transparam <- param_trans(t(gamma_gamma0$gamma) %*% init$beta,
                            emulator::quad.form( init$Sigma,gamma_gamma0$gamma),
                            emulator::quad.form( init$Sigma,gamma_gamma0$gamma0),
                            initA)

  trueA <- as.matrix(find_A_from_gamma(trueparam$Gamma),col=u)
  gamma_gamma0 <- find_gammas_from_A(trueA)

  transparam_true <- param_trans(t(gamma_gamma0$gamma) %*% trueparam$beta,
                            emulator::quad.form( trueparam$Sigma,gamma_gamma0$gamma),
                            emulator::quad.form( trueparam$Sigma,gamma_gamma0$gamma0),
                            trueA)

  Omegatilde <- transparam_true$Omegatilde
  Omega0tilde <- transparam_true$Omega0tilde
  etatilde <- transparam$etatilde
  mu <- transparam_true$mu
  A <- initA
  xi <- rep(1,r-u)

  Yc <- scale(Y,center=TRUE,scale=FALSE)
  Xc <- scale(X,center=TRUE,scale=FALSE)

  rw_var <- rep(1,dim(A)[1])

  postlist <- list()
  for(iter in 1:mcmc.num){
    #A <- initA
    gamma_gamma0 <- find_gammas_from_A(A)
    C1 <- emulator::quad.form(crossprod(Yc) +emulator::quad.form(M,t(B0)) -
                                emulator::quad.form.inv(crossprod(Xc) + M,crossprod(Xc,Yc) + M %*% t(B0)),
                              gamma_gamma0$CA)
    Omegatilde <- CholWishart::rInvWishart(1,df=n+nu+1,Sigma = C1 + psi*gamma_gamma0$CAtCA)[,,1]
    #Omega0tilde <- CholWishart::rInvWishart(1,df=n+nu0+1,
    #                                        Sigma= emulator::quad.form(crossprod(Yc),gamma_gamma0$DA) +
    #                                          psi0*gamma_gamma0$DAtDA)[,,1]

    etabar <- t(solve(crossprod(Xc) + M,crossprod(Xc,Yc) + M %*% t(B0)) %*% gamma_gamma0$CA)
    etatilde <- LaplacesDemon::rmatrixnorm(etabar,LaplacesDemon::as.positive.definite(as.matrix(Omegatilde)),
                                           LaplacesDemon::as.positive.definite(solve(crossprod(Xc) + M)))


    #    muC <- mvnfast::rmvn(1,mu=t(gamma_gamma0$CA) %*% colMeans(Y) - etatilde %*% colMeans(X),
    #                         sigma= Omegatilde/n)
    #    muD <- mvnfast::rmvn(1,mu=t(gamma_gamma0$DA) %*% colMeans(Y),
    #                         sigma= Omega0tilde/n)
    #    mu <- solve(rbind(t(gamma_gamma0$CA),t(gamma_gamma0$DA)),c(muC,muD))



    sqrtCAtCA_inv <- sqrtmatinv(gamma_gamma0$CAtCA)
    sqrtDAtDA_inv <- sqrtmatinv(gamma_gamma0$DAtDA)
    eta <- sqrtCAtCA_inv %*% etatilde
    Omega <- sqrtCAtCA_inv %*% Omegatilde %*% sqrtCAtCA_inv
    Omega0 <- sqrtDAtDA_inv %*% Omega0tilde %*% sqrtDAtDA_inv


    #Sigma <- gamma_gamma0$gamma %*% Omega %*% t(gamma_gamma0$gamma) +
    #  gamma_gamma0$gamma0 %*% Omega0 %*% t(gamma_gamma0$gamma0)
    #mu <- mvnfast::rmvn(1,mu=(barY- gamma_gamma0$gamma %*% eta %*% barX) ,
    #                     sigma= Sigma/n)


    tausqvec <- xi *tau1^2 + (1-xi) * tau0^2

    #gen A
    lpd_val <- lpd_A_pred(A,Y,X,eta,Omega0,Omega,tausqvec,
                          Betashapes,groupindex)
    MHsample <- rwmh_rowwise(lpd_val, rw_var,A,Y,X,eta,
                             Omega0,Omega,tausqvec,
                             Betashapes,groupindex)

    #tune_int <- iter-as.integer(mcmc.num/4)
    tune_int <- iter
    if(tune_int>0) rw_var <- rw_var * exp(tune_int^(-0.7)*(MHsample$alphas-0.44))
    A <- MHsample$A


    for(j in 1:(r-u)){
      logpropratio <- u*log(tau0/tau1) + LaplacesDemon::logit(max(min(0.9999,phi[groupindex[j]]),0.0001))+
        sum(A[j,]^2)/(2*tau0^2)*(1- (tau0/tau1)^2)
      propxi <- LaplacesDemon::invlogit(logpropratio)
      xi[j] <- rbinom(1,1,propxi)

    }

    for(k in 1:K){
      phi[k] <- rbeta(1,Betashapes[k,1]+sum(xi[groupindex==k]),
                      Betashapes[k,2]+sum(1-xi[groupindex==k]))
    }
    postlist[[iter]] <- list(mu=mu,eta=eta,Omega=Omega,
                             Omega0=Omega0,A=A,xi=xi,phi=phi,
                             alphas=MHsample$alphas)

  }


  return(list(postlist=postlist,BIC_MCMC=BIC_MCMC(postlist,X,Y),
              postmean = summary_postsamples_MH(postlist[-(1:as.integer(mcmc.num/2))])))
}


summary_postsamples_MH <- function(mcmc.list){

  xi.mean <- mcmc.list %>% purrr::map(~.x$xi) %>% do.call("rbind",.) %>% colMeans


  beta.mean <- mcmc.list %>%
    purrr::map(~find_gammas_from_A(.x$A)$gamma %*% .x$eta) %>%
    abind::abind(along=3) %>% apply(c(1,2),mean)
  betaxi = mcmc.list %>%
    purrr::map(~find_gammas_from_A(.x$A*.x$xi)$gamma %*% .x$eta) %>%
    abind::abind(along=3) %>% apply(c(1,2),mean)
  list(xi=xi.mean,beta=beta.mean,betaxi=betaxi
  )

}


