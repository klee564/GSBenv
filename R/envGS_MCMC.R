


#'
#' @examples
#' r <- 5
#' p <- 4
#' u <- 2
#' xi.tru <- c(1,1,0)
#' n <- 100
#'
#'
#' groupind <- c(1,1,2)
#' K <- max(groupind)
#' Betashapes <- matrix(1/1000,K,2)
#' group_prior <- list(groupind=groupind,Betashapes=Betashapes)
#'
#' param <- generate_sparsepar(r,p,u, xi.tru)
#'
#' inputdata <- generate_data(param,n)
#' X=inputdata$X;Y=inputdata$Y
#'
#' hyperlist <- list(M = diag(10^(-6),p),B0 = matrix(0,r,p),Psi = diag(1,u),
#' Psi0 = diag(1,r-u),nu = u +1,nu0 = r-u +1,kappa = 10^(-6))
#'
#' Alist <- envGS_MHA(X=inputdata$X,Y=inputdata$Y,u,group_prior,tau0=0.05,tau1=1,hyperlist=hyperlist)
#' Alist[[20]]$A <- param$A.tru
#' postlist <- generate_others(Alist,X,Y,hyperlist)
#'
#'
generate_others <- function(Alist,X,Y,hyperlist){

  n <- nrow(X)
  barY <- colMeans(Y)
  barX <- colMeans(X)

  Xc <- scale(X,scale=FALSE)
  XctXc <- crossprod(Xc)
  barXbarXt <- crossprod(t(barX))

  kappa <- hyperlist$kappa
  M <- hyperlist$M
  Psi <- hyperlist$Psi
  Psi0 <- hyperlist$Psi0
  nu <- hyperlist$nu
  nu0 <- hyperlist$nu0

  SS <- SS_hyper(X,Y,kappa,M)
  Omegafactor <- SS$Omegafactor
  Omega0factor <- SS$Omega0factor
  CXY <- SS$CXY
  etaM <- LaplacesDemon::as.positive.definite(solve(XctXc + n*kappa/(n+kappa)*barXbarXt +M))

  samplingfun <- function(A){
    gamma_gamma0 <- find_gammas_from_A(A)
    gamma <- gamma_gamma0$gamma
    gamma0 <- gamma_gamma0$gamma0

    ## Omega,Omega0 sample
    Omega <-
      CholWishart::rInvWishart(1,n+nu,Psi+ t(gamma) %*% Omegafactor %*% gamma)[,,1]

    Omega0 <-
      CholWishart::rInvWishart(1,n+nu0,Psi0 + t(gamma0) %*% Omega0factor %*%gamma0)[,,1]

    ## eta sample
    etabar <- t(CXY %*% gamma)
    eta <-LaplacesDemon::rmatrixnorm(etabar,Omega ,etaM)

    ## mu sample
    beta <- gamma %*% eta
    Sigma <- gamma %*% Omega %*% t(gamma) + gamma0 %*% Omega0 %*% t(gamma0)
    mu <- mvnfast::rmvn(1,mu=n/(n+kappa)*(barY- beta %*% barX) ,
                        sigma= Sigma/(n+kappa))
    list(mu=mu,eta=eta,Omega=Omega,Omega0=Omega0,beta=beta,Sigma=Sigma)
  }
  burnind <- 1:as.integer(length(Alist)/2)
  otherlist <- purrr::map(Alist[-burnind],~samplingfun(.x$A))
  postlist <- purrr::map2(Alist[-burnind],otherlist,~c(.x,.y))

  list(postlist=postlist,BIC= BIC_MCMC(postlist,X,Y))
}



generate_othermean <- function(A,X,Y,hyperlist){

  u <- ncol(A)
  n <- nrow(X)
  r <- ncol(Y)
  barY <- colMeans(Y)
  barX <- colMeans(X)

  Xc <- scale(X,scale=FALSE)
  XctXc <- crossprod(Xc)
  barXbarXt <- crossprod(t(barX))

  kappa <- hyperlist$kappa
  M <- hyperlist$M
  Psi <- hyperlist$Psi
  Psi0 <- hyperlist$Psi0
  nu <- hyperlist$nu
  nu0 <- hyperlist$nu0

  SS <- SS_hyper(X,Y,kappa,M)
  Omegafactor <- SS$Omegafactor
  Omega0factor <- SS$Omega0factor
  CXY <- SS$CXY
  etaM <- LaplacesDemon::as.positive.definite(solve(XctXc + n*kappa/(n+kappa)*barXbarXt +M))

  gamma_gamma0 <- find_gammas_from_A(A)
  gamma <- gamma_gamma0$gamma
  gamma0 <- gamma_gamma0$gamma0

  ## Omega,Omega0 sample
  Omega <- (Psi+ t(gamma) %*% Omegafactor %*% gamma)/(n+nu-u-1)

  Omega0 <- (Psi0 + t(gamma0) %*% Omega0factor %*%gamma0)/(n+nu0-(r-u))

  eta <- t(CXY %*% gamma)

  ## mu
  beta <- gamma %*% eta
  Sigma <- gamma %*% Omega %*% t(gamma) + gamma0 %*% Omega0 %*% t(gamma0)
  mu <- n/(n+kappa)*(barY- beta %*% barX)
  list(mu=mu,eta=eta,Omega=Omega,Omega0=Omega0,beta=beta,Sigma=Sigma)
}


#'
#' @examples
#' r <- 10
#' p <- 2
#' u <- 2
#' G <- 2
#' params <- generate_GSpar(r,p,u,G)
#'
#' groupind <- params$groupind
#' K <- max(groupind)
#' Betashapes <- matrix(1/1000,K,2)
#' group_prior <- list(groupind=groupind,Betashapes=Betashapes)
#'
#'
#' inputdata <- generate_data(params$param,100)
#' X=inputdata$X;Y=inputdata$Y
#' hyperlist <- list(M = diag(10^(-6),p),B0 = matrix(0,r,p),Psi = diag(1,u),
#' Psi0 = diag(1,r-u),nu = u +1,nu0 = r-u +1,kappa = 10^(-6))
#' postlist <- envGS_MCMC(X=inputdata$X,Y=inputdata$Y,u,group_prior,tau0=0.05,tau1=1,hyperlist=hyperlist)
#'
envGS_MCMC <- function(X,Y,u,group_prior=group_prior,tau0=0.01,tau1=1,hyperlist,mcmc.num=20){
  mcmc.list <- envGS_MHA(X=X,Y=Y,u,group_prior=group_prior,tau0=tau0,tau1=tau1,hyperlist=hyperlist,mcmc.num=mcmc.num) %>%
    generate_others(X=X,Y=Y,hyperlist=hyperlist)

  xi.mean <- mcmc.list$postlist %>% purrr::map(~.x$xi) %>% do.call("rbind",.) %>% colMeans
  beta.mean <- mcmc.list$postlist %>%
    purrr::map(~.x$beta) %>%
    abind::abind(along=3) %>% apply(c(1,2),mean)
  A.mean <- mcmc.list$postlist %>%
    purrr::map(~.x$A) %>%
    abind::abind(along=3) %>% apply(c(1,2),mean)

  #Sigma.mean <- mcmc.list$postlist %>%
  #  purrr::map(~.x$Sigma) %>%
  #  abind::abind(along=3) %>% apply(c(1,2),mean)
  mu.mean <- mcmc.list$postlist %>%
    purrr::map(~.x$mu) %>% do.call("rbind",.) %>% colMeans

  #postximean <- generate_othermean(A.mean * round(xi.mean),X,Y,hyperlist)
  betaxi = mcmc.list$postlist %>%
    purrr::map(~.x$beta * c(rep(1,u),.x$xi)) %>%
    abind::abind(along=3) %>% apply(c(1,2),mean)
  list(A=A.mean,xi=xi.mean,beta=beta.mean,betaxi=betaxi,BIC=mcmc.list$BIC,mu=mu.mean)
       #Sigma = Sigma.mean,mu=mu.mean)
}


