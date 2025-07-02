

#'
#' @examples
#' r <- 10
#' p <- 2
#' u <- 2
#' G <- 2
#' params <- generate_GSpar(r,p,u,G)
#' ex_groupind=c(rep(G+1,u),params$groupind)
#' inputdata <- generate_data(params$param,100)
#'
#' validdata <- generate_data(params$param,20)
#' newX=validdata$X;newY=validdata$Y
#'
#' K <- max(ex_groupind)
#' ex_group_prior <- list(ex_groupind=ex_groupind,
#' ex_Betashapes = matrix(1/1000,K,2))
#' X=inputdata$X;Y=inputdata$Y
#' ex_hyperlist <- list(M = diag(10^(-6),p),B0 = matrix(0,r,p),Psi = diag(1,r),
#' Psi0 = diag(1,r),kappa = 10^(-6))
#' res <- wrap_envGS_MCMC(X,Y,u+1,ex_group_prior=ex_group_prior,
#' tau0=0.05,tau1=10,ex_hyperlist=ex_hyperlist,mcmc.num=20)
#' testerr(res$envres,newX,newY)
#'
testerr <- function(resGS,newX,newY){

  mean((newY - rep(1,nrow(newX)) %*% t(resGS$mu) - newX %*% t(resGS$beta))^2)
}



#'
#' @examples
#' r <- 10
#' p <- 2
#' u <- 2
#' G <- 2
#' params <- generate_GSpar(r,p,u,G)
#' ex_groupind=c(rep(G+1,u),params$groupind)
#' inputdata <- generate_data(params$param,100)
#'
#' validdata <- generate_data(params$param,20)
#' newX=validdata$X;newY=validdata$Y
#'
#' K <- max(ex_groupind)
#' ex_group_prior <- list(ex_groupind=ex_groupind,
#' ex_Betashapes = matrix(1/1000,K,2))
#' X=inputdata$X;Y=inputdata$Y
#' ex_hyperlist <- list(M = diag(10^(-6),p),B0 = matrix(0,r,p),Psi = diag(1,r),
#' Psi0 = diag(1,r),kappa = 10^(-6))
#' tres1 <- CVfun(X,Y,1,tau0=0.05,tau1=10,ex_group_prior,ex_hyperlist,mcmc.num=100)
#' tres2 <- CVfun(X,Y,2,tau0=0.05,tau1=10,ex_group_prior,ex_hyperlist,mcmc.num=100)
#' tres3 <- CVfun(X,Y,3,tau0=0.05,tau1=10,ex_group_prior,ex_hyperlist,mcmc.num=100)
#'
CVfun <- function(X,Y,u,tau0,tau1,ex_group_prior,ex_hyperlist,mcmc.num){

  n <- nrow(X)
  validx <- caret::createMultiFolds(1:n,k=10,times=5)

  library(furrr)
  plan(multisession, workers = 10)


  trainfit <- validx %>%
    furrr::future_map(~wrap_envGS_MCMC(X[.x,],Y[.x,],u,ex_group_prior=ex_group_prior,
                                tau0=tau0,tau1=tau1,ex_hyperlist=ex_hyperlist,mcmc.num=mcmc.num))


  purrr::map(1:length(trainfit),
              ~testerr(trainfit[[.x]]$envres,newX=X[-validx[[.x]],],newY=Y[-validx[[.x]],])) %>%
    unlist %>% mean

}

#'
#' @examples
#' r <- 10
#' p <- 2
#' u <- 2
#' G <- 2
#' params <- generate_GSpar(r,p,u,G)
#' ex_groupind=c(rep(G+1,u),params$groupind)
#' inputdata <- generate_data(params$param,100)
#'
#' K <- max(ex_groupind)
#' ex_group_prior <- list(ex_groupind=ex_groupind,
#' ex_Betashapes = matrix(1/1000,K,2))
#' X=inputdata$X;Y=inputdata$Y
#' ex_hyperlist <- list(M = diag(10^(-6),p),B0 = matrix(0,r,p),Psi = diag(1,r),
#' Psi0 = diag(1,r),kappa = 10^(-6))
#' tres1 <- wrap_CV(X,Y,uvec=1:3,tau0vec=c(0.01,0.05,0.1),tau1=1,ex_group_prior,ex_hyperlist,mcmc.num=100)
#' zz <- wrap_envGS_MCMC(X,Y,2,ex_group_prior=ex_group_prior,tau0=0.05,tau1=10,ex_hyperlist=ex_hyperlist,mcmc.num=200)
#'
#'
wrap_CV <- function(X,Y,uvec,tau0vec,tau1,ex_group_prior,ex_hyperlist,mcmc.num){
  library(magrittr)
  tdf <- expand.grid(u=uvec,tau0=tau0vec)

  value <- 1:nrow(tdf) %>%
      purrr::map(~CVfun(X=X,Y=Y,u=tdf$u[.x],tau0=tdf$tau0[.x],tau1=tau1,ex_group_prior=ex_group_prior,
                        ex_hyperlist=ex_hyperlist,mcmc.num=mcmc.num)) %>% unlist

  data.frame(tdf,value=value)

}







