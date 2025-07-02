
#'
#' @examples
#' r <- 5
#' p <- 4
#' u <- 2
#' xi.tru <- c(1,1,0)
#'
#' groupind <- c(1,1,2)
#' K <- max(groupind)
#' Betashapes <- matrix(1/1000,K,2)
#' group_prior <- list(groupind=groupind,Betashapes=Betashapes)
#'
#' param <- generate_sparsepar(r,p,u, xi.tru)
#'
#' inputdata <- generate_data(param,50)
#' X=inputdata$X;Y=inputdata$Y
#' hyperlist <- list(M = diag(10^(-6),p),B0 = matrix(0,r,p),Psi = diag(1,u),
#' Psi0 = diag(1,r-u),nu = u +1,nu0 = r-u +1,kappa = 10^(-6))
#' postlist <- envGS_MHA(X=inputdata$X,Y=inputdata$Y,u,group_prior,tau0=0.01,tau1=100,hyperlist=hyperlist)
#'
envGS_MHA <- function(X,Y,u,group_prior=group_prior,tau0,tau1,hyperlist,mcmc.num=200){

  postlist <- list()

  r <- ncol(Y)
  p <- ncol(X)
  n <- nrow(X)

  groupindex <- group_prior$groupind
  K <- max(groupindex)
  phi <- rep(1/2,K)

  #set hyperparameter
  M <- hyperlist$M
  B0 <- hyperlist$B0
  Psi <- hyperlist$Psi
  Psi0 <- hyperlist$Psi0
  nu <- hyperlist$nu
  nu0 <- hyperlist$nu0
  kappa <- hyperlist$kappa

  Betashapes <- group_prior$Betashapes

  #set init
  init <- Renvlp::env(X=X,Y=Y,u = u)
  initA <- as.matrix(find_A_from_gamma(init$Gamma),col=u)
  gamma_gamma0 <- find_gammas_from_A(initA)

  A <- initA
  xi <- rep(1,r-u)
  rw_var <- rep(1,dim(A)[1])

  SS <- SS_hyper(X,Y,kappa,M)

  postlist <- list()
  for(iter in 1:mcmc.num){
    tausqvec <- xi *tau1^2 + (1-xi) * tau0^2

    #gen A
    lpd_val <- lpd_A_pred(n,A,SS$Omegafactor,SS$Omega0factor,Psi,Psi0,nu,nu0,tausqvec)
    MHsample <- rwmh_rowwise(lpd_val, rw_var,n,A,SS$Omegafactor,SS$Omega0factor,Psi,Psi0,nu,nu0,tausqvec)

    #tune_int <- iter-as.integer(mcmc.num/4)
    tune_int <- iter
    if(tune_int>0) rw_var <- rw_var * exp(tune_int^(-0.7)*(MHsample$alphas-0.44))
    A <- MHsample$A


    for(k in 1:K){

      # if (autotune_size == "single") tau_curr <- tau
      xig_star <- rbinom(sum(groupindex==k),1,1/2)
      xig_pre <- xi[groupindex==k]


      tausqvec_star <- xig_star *tau1^2 + (1-xig_star) * tau0^2
      tausqvec_pre <- xig_pre *tau1^2 + (1-xig_pre) * tau0^2

      lpd_xi_star <- 1:length(tausqvec_star) %>%
        purrr::map(~mvnfast::dmvn(A[which(groupindex==k)[.x],],mu=rep(0,u),
                                  sigma = diag(tausqvec_star[.x],u),log = T)) %>% unlist %>%sum

      lpd_xi_pre <- 1:length(tausqvec_pre) %>%
        purrr::map(~mvnfast::dmvn(A[which(groupindex==k)[.x],],mu=rep(0,u),
                                  sigma = diag(tausqvec_pre[.x],u),log = T)) %>% unlist %>% sum

      lpd_xi_star <- lpd_xi_star +
        extraDistr::dbbinom(sum(xig_star),length(xig_star),
                            Betashapes[k,1],Betashapes[k,2],log = TRUE)
                            #sum(xig_star)+ Betashapes[k,1],sum(1-xig_star)+Betashapes[k,2])
      lpd_xi_pre <- lpd_xi_pre +
        extraDistr::dbbinom(sum(xig_pre),length(xig_pre),
                            Betashapes[k,1],Betashapes[k,2],log = TRUE)
                            #sum(xig_pre)+Betashapes[k,1],sum(1-xig_pre)+Betashapes[k,2])

      aprob <- min(exp(lpd_xi_star-lpd_xi_pre),1)
      if (runif(1) < aprob) {
        # accept
        xi[groupindex==k] <- xig_star
      }
    }
    #for(j in 1:(r-u)){
    #  #logpropratio <- u*log(tau0/tau1) + LaplacesDemon::logit(max(min(0.999,phi[groupindex[j]]),0.0001))+
    #  logpropratio <- u*log(tau0/tau1) + max(min(LaplacesDemon::logit(phi[groupindex[j]]),10^(300)),-10^(300))+
    #    sum(A[j,]^2)/(2*tau0^2)*(1- (tau0/tau1)^2)
    #  propxi <- LaplacesDemon::invlogit(logpropratio)
    #  xi[j] <- rbinom(1,1,propxi)
    #}
    for(k in 1:K){
      phi[k] <- rbeta(1,Betashapes[k,1]+sum(xi[groupindex==k]),
                      Betashapes[k,2]+sum(1-xi[groupindex==k]))
    }
    postlist[[iter]] <- list(A=A,xi=xi,phi=phi,alphas=MHsample$alphas)
  }
  return(postlist)
}

