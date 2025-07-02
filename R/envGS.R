
#'
#' @examples
#' r <- 5
#' p <- 4
#' u <- 2
#' xi.tru <- c(1,1,0)
#' groupind <- c(1,1,2)
#'
#' param <- generate_sparsepar(r,p,u, xi.tru)
#'
#' inputdata <- generate_data(param,200)
#' postlist <- envGS(X=inputdata$X,Y=inputdata$Y,u,groupind,mcmc.num=100,hyper=NULL)
#' res <- purrr::map(postlist,~param_invtrans(.x$mu,.x$etatilde,.x$Omegatilde,.x$Omega0tilde,.x$A))
#' mean(((res %>% purrr::map(~.x$beta) %>% purrr::reduce(`+`))/100- param$beta.tru)^2)
#' mean((Renvlp::env(X=inputdata$X,Y=inputdata$Y,u = u)$beta- param$beta.tru)^2)
#'
#' mean(((res %>% purrr::map(~.x$mu) %>% purrr::reduce(`+`))/100- param$mu.tru)^2)
#' mean((Renvlp::env(X=inputdata$X,Y=inputdata$Y,u = u)$mu- param$mu.tru)^2)
#'
envGS <- function(X,Y,u,groupind,mcmc.num=20,hyper=NULL){

  postlist <- list()

  r <- ncol(Y)
  p <- ncol(X)
  n <- nrow(X)

  groupindex <- as.integer(as.factor(groupind))
  K <- max(groupindex)
  phi <- rep(1/2,K)


  #set hyperparameter
  M <- diag(10^(-6),p)
  Msqrt <- sqrtmat(M)
  B0 <- matrix(0,r,p)
  psi <- 1
  psi0 <- 1
  nu <- u +1
  nu0 <- r-u +1
  tau0 <- sqrt(r/n)
  tau1 <- 1

  Betashapes <- matrix(1/1000,K,2)


  #set init
  init <- Renvlp::env(X=X,Y=Y,u = u)
  initA <- find_A_from_gamma(init$Gamma)
  gamma_gamma0 <- find_gammas_from_A(initA)
  transparam <- param_trans(t(gamma_gamma0$gamma) %*% init$beta,
                            emulator::quad.form( init$Sigma,gamma_gamma0$gamma),
                            emulator::quad.form( init$Sigma,gamma_gamma0$gamma0),
                            initA)
  Omegatilde <- transparam$Omegatilde
  Omega0tilde <- transparam$Omega0tilde
  etatilde <- transparam$etatilde
  mu <- init$mu
  A <- initA
  xi <- rep(1,r-u)

  Yc <- scale(Y,center=TRUE,scale=FALSE)
  Xc <- scale(X,center=TRUE,scale=FALSE)


  postlist <- list()
  for(iter in 1:mcmc.num){
    #A <- initA
    gamma_gamma0 <- find_gammas_from_A(A)
    C1 <- emulator::quad.form(crossprod(Yc) +emulator::quad.form(M,t(B0)) -
                                emulator::quad.form.inv(crossprod(Xc) + M,crossprod(Xc,Yc) + M %*% t(B0)),
                              gamma_gamma0$CA)
    Omegatilde <- CholWishart::rInvWishart(1,df=n+nu+1,Sigma = C1 + psi*gamma_gamma0$CAtCA)[,,1]
    Omega0tilde <- CholWishart::rInvWishart(1,df=n+nu0+1,
                                            Sigma= emulator::quad.form(crossprod(Yc),gamma_gamma0$DA) +
                                              psi0*gamma_gamma0$DAtDA)[,,1]

    etabar <- t(solve(crossprod(Xc) + M,crossprod(Xc,Yc) + M %*% t(B0)) %*% gamma_gamma0$CA)
    etatilde <- LaplacesDemon::rmatrixnorm(etabar,LaplacesDemon::as.positive.definite(Omegatilde),
                                           LaplacesDemon::as.positive.definite(solve(crossprod(Xc) + M)))

    #muC <- mvnfast::rmvn(1,mu=t(gamma_gamma0$CA) %*% colMeans(Y) - etatilde %*% colMeans(X),
    #                     sigma= Omegatilde/n)
    #muD <- mvnfast::rmvn(1,mu=t(gamma_gamma0$DA) %*% colMeans(Y),
    #                     sigma= Omega0tilde/n)
    #mu <- solve(rbind(t(gamma_gamma0$CA),t(gamma_gamma0$DA)),c(muC,muD))

    #mu <- rep(0,r)
    #etatilde <- transparam$etatilde
    #Omegatilde <- transparam$Omegatilde
    #Omega0tilde <- transparam$Omega0tilde


    #gen A
    Ymu <- sweep(Y,2,mu)
    Omegatildeinv <- solve(Omegatilde)
    Omega0tildeinv <- solve(Omega0tilde)

    for(j in 1:(r-u)){
      tausq <- xi[j] * tau1^2 + (1-xi[j]) * tau0^2
      for(k in 1:u){
        A[j,k] <- update_a(j,k, n+as.integer((nu+nu0)/2),
                           Ymu,X,Omegatildeinv,Omega0tildeinv,etatilde,
                           A=A,tausq=tausq,B0=B0,Msqrt=Msqrt,psi=psi,psi0=psi0,Miter=100)
      }
      logpropratio <- u*log(tau0/tau1) + LaplacesDemon::logit(phi[groupindex[j]])+
        sum(A[j,]^2)/(2*tau0^2)*(1- (tau0/tau1)^2)
      propxi <- LaplacesDemon::invlogit(logpropratio)
      xi[j] <- rbinom(1,1,propxi)
    }
    for(k in 1:K){
      phi[k] <- rbeta(1,Betashapes[k,1]+sum(xi[groupindex==k]),
                      Betashapes[k,2]+sum(1-xi[groupindex==k]))
    }
    postlist[[iter]] <- list(mu=mu,etatilde=etatilde,Omegatilde=Omegatilde,
                             Omega0tilde=Omega0tilde,A=A,xi=xi,phi=phi)

  }

  return(postlist)
}


#'
#' @examples
#'
#' r <- 6
#' p <- 4
#' u <- 2
#' xi.tru <- sample(c(0,1),r-u,replace=TRUE)
#'
#' param <- generate_sparsepar(r,p,u, xi.tru)
#'
#' inputdata <- generate_data(param,100)
#' mcmc.list <- envGS(X=inputdata$X,Y=inputdata$Y,u,mcmc.num=100,hyper=NULL)
#'
summary_postsamples <- function(mcmc.list){

  xi.mean <- mcmc.list %>% purrr::map(~.x$xi) %>% do.call("rbind",.) %>% colMeans


  beta.mean <- mcmc.list %>% purrr::map(~param_trans_inv(.x)$beta) %>%
    abind::abind(along=3) %>% apply(c(1,2),mean)

  list(xi=xi.mean,beta=beta.mean)

}



