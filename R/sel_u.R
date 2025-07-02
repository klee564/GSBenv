
#'
#' @examples
#' r <- 10
#' p <- 2
#' u <- 2
#' G <- 2
#' nvec <- c(100)
#' tau0vec=c(0.01,0.05,0.1)
#' uvec <- 1:3
#' params <- generate_GSpar(r,p,u,G)
#' inputdata <- generate_data(params$param,50)
#'
#' G <- max(params$groupind)
#' ex_Betashapes <- matrix(c(0.00001,0.00001),G+1,2,byrow=TRUE)
#' ex_group_prior <- list(ex_groupind=c(rep(G+1,u),params$groupind),ex_Betashapes=ex_Betashapes)
#' X=inputdata$X;Y=inputdata$Y
#' uvec <- c(1,2,3)
#' tau0vec <- c(0.01,0.05)
#' ex_hyperlist <- list(M = diag(10^(-6),p),B0 = matrix(0,r,p),Psi = diag(1,r),
#' Psi0 = diag(1,r),kappa = 10^(-6))
#' rr <- envGS_selu(Y,X,ex_group_prior,uvec,tau0vec,tau1=1,ex_hyperlist=ex_hyperlist)
#'
envGS_selu <- function(Y,X,ex_group_prior,uvec,tau0vec,tau1=1,ex_hyperlist=ex_hyperlist, ncore=1){
  library(magrittr)
  tdf <- expand.grid(u=uvec,tau0=tau0vec)

  if(ncore>1){
    library(future)
    plan(multisession, workers = ncore)
    res <- 1:nrow(tdf) %>%
      furrr::future_map(~wrap_envGS_MCMC(Y=Y,X=X,u=tdf$u[.x],
                                       ex_group_prior=ex_group_prior,
                                       tau0 = tdf$tau0[.x],tau1=tau1,ex_hyperlist=ex_hyperlist,mcmc.num=1000))
  } else{
    res <- 1:nrow(tdf) %>%
      purrr::map(~wrap_envGS_MCMC(Y=Y,X=X,u=tdf$u[.x],
                                       ex_group_prior=ex_group_prior,
                                       tau0 = tdf$tau0[.x],tau1=tau1,ex_hyperlist=ex_hyperlist,mcmc.num=1000))
  }

  optind <- which.min(res %>% purrr::map(~.x$envres$BIC) %>% unlist)

  list(result=res[[optind]],u= tdf$u[optind],tau0= tdf$tau0[optind],
       BIC=data.frame(tdf,BIC=res %>% purrr::map(~.x$envres$BIC) %>% unlist))
}


#'
#' @examples
#' r <- 10
#' p <- 2
#' u <- 2
#' G <- 2
#' uvec <- 1:3
#' params <- generate_GSpar(r,p,u,G)
#' inputdata <- generate_data(params$param,50)
#' X <- inputdata$X
#' Y <- inputdata$Y
#' zz <- senv_selu(Y,X,uvec = 2,lambda1=0,lambda2=1)
#'
senv_selu <- function(Y,X,uvec, ncore=1,init=NULL,lambda1=NULL,lambda2=NULL){

  library(magrittr)
  tdf <- expand.grid(u=uvec)
  if(ncore>1){
    library(future)
    plan(multisession, workers = ncore)
    res <- 1:nrow(tdf) %>%
      furrr::future_map(~spenv(X=X,Y=Y,u=tdf$u[.x],init=init,lambda1=lambda1,lambda2=lambda2))
  } else{
    res <- 1:nrow(tdf) %>%
      purrr::map(~spenv(X=X,Y=Y,u=tdf$u[.x],init=init,lambda1=lambda1,lambda2=lambda2))
  }


  optind <- which.min(res %>% purrr::map(~.x$BIC) %>% unlist)
  list(result=res[[optind]],u= tdf$u[optind],
       BIC=data.frame(tdf,BIC=res %>% purrr::map(~.x$BIC) %>% unlist))
}
